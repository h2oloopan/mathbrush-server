#ifdef _DEBUG
	#undef _DEBUG
	#include <Python.h>
	#define _DEBUG
#else
	#include <Python.h>
#endif

#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <errno.h>
#include <fcntl.h>
#include <cstdlib>
#include <cstdio>
#include <signal.h>
#include <cstring>
#include <cassert>

#include <set>
#include <string>
#include <string.h>
#include <sstream>
#include <stdexcept>

#include "cmdcode.h"
#include "expr.h"
#include "io.h"
#include "log.h"
#include "lconv.h"
#include "MathRecognizer.h"

static pid_t childpid = -1;

// Directory containing the Python files
static const std::string CASERVER_DIR = "/home/rmprosse/mathBrushResearchII/Code/C++/mathreco/TRUNK/src/caserver/"; 

static void
killreco() {
	//clear_relations();
	scg::ShutdownRecognizer();
}

static void
hup() {
	killreco();
	int rep = CASREP_HUP;
	write(1, (void *)&rep, sizeof(rep));
	exit(-1);
	//flushandexit(rxfd, -1);
}

#define ESTACK_MAX 32
static exprtree *estack[ESTACK_MAX];
static size_t estacki = 0;
#define ESTACK_SIZE() estacki

static exprtree *
currexpr(void) {
	if (estacki == 0) {
		return 0;
	}
	return estack[estacki-1];
}

static int
pushexpr(exprtree *expr) {
	if (estacki == ESTACK_MAX) {
		return CASREP_TOOLARGE;
	}
	else {
		estack[estacki] = expr;
		++estacki;
		return CASREP_OK;
	}
}

static exprtree *
popexpr() {
	if (estacki == 0) {
		return 0;
	}
	else {
		--estacki;
		return estack[estacki];
	}
}

static int
cmdexpr(void) {
	int e;
	exprtree *tree;
	e = readexpr(0, &tree);
	if (e == CASREP_OK) {
		pushexpr(tree);
		logtime();
		std::string s;
		tree->writesage(s);
		logmsg("Loaded expr %s\n", s.c_str());
	}

	return e;
}

static std::set<std::string> declnames;

// Name of Python variable which stores the returned expression 
// (arbitrary, but shared between writesagecmd and invokepythonsage)
static const std::string EXPR_VAR_NAME = "expr"; 

static PyObject *pMainModule; // Python main module
static PyObject *pExprTreeClass; // Python ExpressionTree class

/* 
 * Writes the code for a specified Sage command to a string. The Sage code is determined by:
 * @cmdlead which specifies the leading code (e.g. "factor("),
 * @cmdtrail which specifies the trailing code (e.g. ")"), and
 * @nargs which specifies the number of arguments to pass to the Sage function.
 * The arguments are expression trees which are popped from the stack.
 * The resulting Sage code is written into the reference @code.
 * Returns the reply code (defined in cmdcode.h).
 */
static int
writesagecmd(std::string& code, const char *cmdlead, const char *cmdtrail, size_t nargs) {
	
	code = "";
	
	#define MAXARGS 4
	static exprtree *args[MAXARGS];
	int e = CASREP_OK;
	size_t i;
	exprtree *tree;
	if (nargs > MAXARGS) {
		return CASREP_SYSERR;
	}
	if (ESTACK_SIZE() < nargs + 1) {
		return CASREP_MSGERR;
	}

	for (i = nargs; i > 0; --i) {
		args[i-1] = popexpr();
		assert(args[i-1]);
	}
	tree = popexpr();
	assert(tree);
	
	std::string vardecls = "";
	declnames.clear();
	e = tree->writevardecls(vardecls, declnames);
	if (e != CASREP_OK) {
		// may have written partial data to sage, so provoke a keyboard interrupt
		kill(childpid, SIGINT);
 		logmsg("writing tree failed\n");
		return e;
	}
	
	code += vardecls + "\n\n";

	std::string sagecode = std::string(cmdlead);
	e = tree->writesage(sagecode);
	if (e == CASREP_OK) {
		for (i = 0; i < nargs; ++i) {
			sagecode += ", ";
			e = args[i]->writesage(sagecode);
			if (e != CASREP_OK) {
				break;
			}
		}
	}
	if (e == CASREP_OK)
		sagecode += std::string(cmdtrail);

	logmsg("Generated Sage code: %s\n", sagecode.c_str());

	// invoke Sage's preparser (see: ~/Sage/sage-4.8/devel/sage/sage/misc/preparser.py)
	code += "try: " + EXPR_VAR_NAME + " = ExpressionTree(eval(preparse(\"" + sagecode + "\")))\n";
	code += "except (ValueError): " + EXPR_VAR_NAME + " = \"BAD\"\n";
	code += "except (TypeError, AttributeError, ArithmeticError): " + EXPR_VAR_NAME + " = \"CMD\"\n";
	code += "except: " + EXPR_VAR_NAME + " = \"SYS\"\n";
	
	return e;
	
} // writesagecmd

/*
 * A helper class for "guarding" references for the Python C API.
 */
class Guard {
	PyObject *obj;
public:
	Guard(PyObject *obj_): obj(obj_) {
		if (!obj)
			throw std::runtime_error("request to guard NULL python object\n");
	}
	~Guard() {
		Py_DECREF(obj);
	}
};

/* 
 * Converts an expression tree from Python to C++.
 * @py_tree represents a Python instance of ExpressionTree.
 * Returns a pointer to a corresponding instance of exprtree.
 */
static exprtree* 
convertexprtree(PyObject* py_tree)
throw (std::runtime_error) {
	
	exprgrp* cpp_tree;
	long sid;
	
	PyObject* sid_attr = PyString_FromString("sid");
	Guard sid_attr_guard(sid_attr);
	if (PyObject_HasAttr(py_tree, sid_attr)) { // todo: does this check actually work? (always true)
		
		PyObject* id = PyObject_GetAttr(py_tree,sid_attr);
		Guard id_guard(id);      
		if (id && PyInt_Check(id)) { 
			
			sid = PyInt_AsLong(id);
			if (sid == scg::TERMINAL_EXPR) {
				
				PyObject* terms_attr = PyString_FromString("terms");
				Guard terms_attr_guard(terms_attr);
				if (PyObject_HasAttr(py_tree, terms_attr)) {
					
					PyObject* terms = PyObject_GetAttr(py_tree,terms_attr);
					Guard terms_guard(terms);      
					if (terms && PyString_Check(terms)) { 
						
						exprleaf* leaf = new exprleaf();
						const char* t = PyString_AsString(terms);
						
						// check for special symbols by name (e.g. pi, leq)
						static const scg::symbol *s;
						if (s = scg::symdb_findsymbol_name(std::string(t))) 
						{ 
							leaf->addterm(s->unicode); // use unicode
						} 
						else
						{
							int size = strlen(t);
							for (int i=0; i<size; i++)
								leaf->addterm(t[i]);
						}

						return leaf; 	  
					} 
					else {
						PyErr_Format(PyExc_TypeError, "unrecognized %s object!", Py_TYPE(terms)->tp_name);
						throw std::runtime_error("unrecognized terminal object in Python expression tree\n");	          
					}
					
				}
			} 
			else  
				cpp_tree = new exprgrp(sid); 
		}
	}
	
	if (cpp_tree == NULL)
		throw std::runtime_error("Unable to get semantic ID from python object\n");
	
	PyObject* list = PyObject_GetAttrString(py_tree, "children");
	if (!list && PyErr_ExceptionMatches(PyExc_AttributeError)) 
	{		
		PyErr_Clear();  // hasattr does this exact check
		return cpp_tree;
	}
	Guard list_guard(list);
	
	// Construct the children of the expression tree
	Py_ssize_t size = PyList_Size(list);
	for (int i=0; i<size; i++) 
	{
		PyObject* py_child = PyList_GetItem(list,i);
		Py_XINCREF(py_child);
		Guard py_child_guard(py_child);

		if (PyObject_IsInstance(py_child, pExprTreeClass)) { //PyClass_Check(py_child)
			exprtree* cpp_child = convertexprtree(py_child);
			cpp_tree->addchild(cpp_child);
		} 
		else {
			PyErr_Format(PyExc_TypeError, "unrecognized %s object!", Py_TYPE(py_child)->tp_name);
			throw std::runtime_error("unrecognized object from python (expected ExpressionTree)\n");  
		}
	}
	 
	return cpp_tree;
	
} // convertexprtree

/*
 * Invoke Sage via the Python C API and return the expression tree.
 * @code is a reference to the Python/Sage code to be invoked.
 * @tree is a reference to the expression tree converted from the obtained Python instance.
 * Returns a reply code (defined in cmdcode.h).
 */
static int 
invokepythonsage(std::string& code, exprtree*& tree)
{
	try {
		int e;
		tree = NULL;
		
		PyObject *pDict = PyModule_GetDict(pMainModule);

		// execute Python code
		PyRun_String(code.c_str(), Py_file_input, pDict, pDict); 
 		
		// get resulting value
		PyObject *pValue = PyRun_String(EXPR_VAR_NAME.c_str(), Py_eval_input, pDict, pDict); 
		Guard value_guard(pValue);
		
		if (pValue != NULL) 
		{
			// if pValue is instance of ExpressionTree (as it should be if no errors)
			if (PyObject_IsInstance(pValue, pExprTreeClass)) {
				
				tree = convertexprtree(pValue);
				e = CASREP_OK;
			}
			// if something went wrong a string may have been returned
			else if (PyString_Check(pValue)) {
				
				char* s = PyString_AsString(pValue);
				if (!strcmp(s, "BAD"))
					e = CASREP_BADEXPR;
				else if (!strcmp(s, "CMD"))
					e = CASREP_INVALIDCMD;
/*				else if (!strcmp(s, "SYS"))
					e = CASREP_SYSERR;*/
				else
					e = CASREP_SYSERR;
			}
			else
				e = CASREP_SYSERR;
		}
		else 
		{
			PyErr_Print();
			e = CASREP_SYSERR;
		}
		
		return e;
	}
	catch (std::exception& ex) {
		logmsg(ex.what());
		return CASREP_SYSERR;
	}  
} // invokepythonsage

/*
 * Performs a specified Sage command by generating the Sage/Python code (see writesagecmd)
 * and invoking it via embedded Python (see invokepythonsage).
 * The Sage code is determined by:
 * @cmdlead which specifies the leading code (e.g. "factor("),
 * @cmdtrail which specifies the trailing code (e.g. ")"), and
 * @nargs which specifies the number of arguments to pass to the Sage function.
 * The arguments are expression trees which are popped from the stack.
 * The resulting Sage code is written into the reference @code.
 * Writes the resulting tree back to the client and returns a reply code (defined in cmdcode.h).
 */
static int
dosagecmd(const char *cmdlead, const char *cmdtrail, size_t nargs) {
	int e;
	exprtree *tree;

	logtime();

	std::string code;
	if ((e = writesagecmd(code, cmdlead, cmdtrail, nargs)) != CASREP_OK)
		return e;

	if ((e = invokepythonsage(code, tree)) != CASREP_OK) 
		return e;

	if (!tree)
		return CASREP_SYSERR;

	
	std::string s;
	tree->writesage(s);
	logmsg("Converted result to expression tree: %s\n", s.c_str());
	
	std::vector<char> buf;
	e = tree->writestream(buf);
	if (e != CASREP_OK)
		return e;
	
	int32_t rep = CASREP_EXPR;
	logmsg("replying with code %d and tree response of %u bytes\n", rep, buf.size());
	if ((e = fullwrite(1, (void *)&rep, sizeof(rep))) != sizeof(rep))
		return e;
	
	uint32_t wrlen32 = buf.size();
	if ((e = fullwrite(1, (void *)&wrlen32, sizeof(wrlen32))) != sizeof(wrlen32))
		return e;
	
	if ((e = fullwrite(1, (void *)&buf[0], buf.size())) != buf.size())
		return e;
	
	return CASREP_OK;
} // dosagecmd

static int
cmdeval() {
	logmsg("evalexpr\n");
	return dosagecmd("", "", 0);
}

static int
cmdevalnum() {
	logmsg("evalnum\n");
//	return dosagecmd(casfp, casrx, "", ".n()", 0);
	return dosagecmd("numerical_approx(", ")", 0);
}

static int
cmdsimplify() {
	logmsg("simplify\n");
	//return dosagecmd("simplify(", ")", 0);
	return dosagecmd("DoSimplify(", ")", 0); // wrapper
}

static int
cmdexpand() {
	logmsg("expand\n");
	return dosagecmd("expand(", ")", 0);
}

static int
cmdsolve() {
	logmsg("solve\n");
	return dosagecmd("DoSolve(", ")", 1); // Python wrapper
}

static int
cmdsolvesystem() {
	logmsg("solvesystem\n");
	
	// get current expression
	exprtree *tree = currexpr();
	assert(tree);
	
	// extract the variables
	std::string cmdtrail = "";
	std::string vardecls = "";
	declnames.clear();
	int e = tree->writevardecls(vardecls, declnames);
	if (e != CASREP_OK) return e;		
	for (std::set<std::string>::const_iterator i = declnames.begin(); i != declnames.end(); ++i)
		cmdtrail += ", " + *i;
	cmdtrail += ")";
	
	return dosagecmd("solve(", cmdtrail.c_str(), 0); // Python wrapper
}

static int
cmdfactor() {
	logmsg("factor\n");
	return dosagecmd("DoFactor(", ")", 0); // Python wrapper
}

static int
cmddeterminant() {
	logmsg("determinant\n");
	return dosagecmd("det(", ")", 0);
}

static int
cmdinverse() {
	logmsg("inverse\n");
	return dosagecmd("", ".inverse()", 0);
}

static int
cmdrank() {
	logmsg("rank\n");
	return dosagecmd("", ".rank()", 0);
}

static int
cmdnullspace() {
	logmsg("nullspace\n");
	return dosagecmd("", ".right_kernel()", 0);
}

static void
runserver() {
	int e;

	log_init();

	sockaddr_in addr;
	socklen_t len = sizeof(addr);
	e = getpeername(0, (sockaddr *)&addr, &len);
	logtime();
	if (e == 0)
		logmsg("starting caserver for client %s\n", inet_ntoa(addr.sin_addr));
	else
		logmsg("starting caserver for unknown client\n");

	// sleep to give you a chance to connect the debugger
/*	int seconds = 30;
	usleep(seconds * pow(10,6));*/
	
	scg::DisableRelationFiltering();
	scg::DisableSubtreeFiltering();
	//scg::SetGrammar("server");
	//reset_relations();
	scg::InitializeRecognizer();

	for (;;) {
		int32_t cmd;
		int32_t rep;
		e = fullread(0, (void *)&cmd, sizeof(cmd));

		logmsg("\n");
		logtime();
		logmsg("fullread returned %d\n", e);
		
		if (e == CASCMD_HUP)
			break;
		else if (e < 0)
			hup();

		logmsg("read command %d\n", cmd);

		switch (cmd) {
		case CASCMD_HUP: 
			goto die;
		case CASCMD_EXPR:
			rep = cmdexpr();
			break;
		case CASCMD_EVAL:
			rep = cmdeval();
			break;
		case CASCMD_EVALNUM:
			rep = cmdevalnum();
			break;
		case CASCMD_SIMPLIFY:
			rep = cmdsimplify();
			break;
		case CASCMD_EXPAND:
			rep = cmdexpand();
			break;
		case CASCMD_SOLVE:
			rep = cmdsolve();
			break;
		case CASCMD_SOLVESYSTEM:
			rep = cmdsolvesystem();
			break;
		case CASCMD_FACTOR:
			rep = cmdfactor();
			break;
		case CASCMD_DETERMINANT:
			rep = cmddeterminant();
			break;
		case CASCMD_INVERSE:
			rep = cmdinverse();
			break;
		case CASCMD_RANK:
			rep = cmdrank();
			break;
		case CASCMD_NULLSPACE:
			rep = cmdnullspace();
			break;
		default:
			logmsg("invalid command code %s\n", cmd);
			rep = CASREP_HUP;
			break;
		}
		logmsg("got reply %d\n", rep);

		switch (rep) {
		case CASREP_HUP:
			fullwrite(1, (void *)&rep, sizeof(rep));
			goto die;
		default:
			if (rep != CASREP_OK || cmd == CASCMD_EXPR) {
				fullwrite(1, (void *)&rep, sizeof(rep));
			}
			break;
		}
	} // for

die:
	logmsg("client disconnected\n");

	killreco();

	/* wait for sage process to end. SIGCHLD handler will exit this process. */
	int status;
	wait(&status);
}

static void 
initpythonsage(int argc, char *argv[])
{
	// Initialize Python C API
	Py_Initialize();    
	PySys_SetArgv(argc, argv); // rmprosse

	// Change to executing directory
	PyRun_SimpleString("import os");	
	std::string cd = "os.chdir(\"" + CASERVER_DIR + "\")"; 
	PyRun_SimpleString(cd.c_str());

	// import the Sage wrapper module (includes Sage libraries)
	PyRun_SimpleString("from SageWrapper import *");
	
	pMainModule = PyImport_AddModule("__main__");
	//PyObject *module = PyImport_ImportModule("SageWrapper"); 
	//Guard module_guard(module);
	
	PyObject *module = PyImport_AddModule("SageWrapper");	
	pExprTreeClass = PyObject_GetAttrString(module, "ExpressionTree"); 
	Guard cls_guard(pExprTreeClass);	
}

static void 
endpythonsage()
{
	// Clean up
	Py_DECREF(pMainModule);
	Py_Finalize();	
}

int
main(int argc, char *argv[]) {

	initpythonsage(argc, argv);	
	runserver();	
	endpythonsage();

	return 0;
}
