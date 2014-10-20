#include <sys/types.h>
#ifdef WIN32
#include <io.h>
#define NOMINMAX
#include <winsock.h>
#else
#include <sys/socket.h>
#include <netdb.h>
#include <unistd.h>
#endif
#include <errno.h>
//#include <fcntl.h>
#include "io.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <set>

#include "int32.h"
#include "caclient.h"
#include "CASOperation.h"
#include "cmdcode.h"
#include "etype.h"
#include "expr.h"
#include "MathRecoTypes.h"
#include "grammar-values.h"
#include "mathrecognizer-private.h"
#include "symbols.h"
#include "annotate.h"
#include "verb.h"

using scg::verb_out;

//#define CASERVER "scg40.cs.uwaterloo.ca"
#define CASERVER "129.97.170.210"
#define CASERVER_PORT 22222

// reco symbols inserted by exprtree_t -> ExpressionTree conversion
static const scg::symbol *S_infin = 0;
static const scg::symbol *S_plus = 0;
static const scg::symbol *S_horzline = 0;
static const scg::symbol *S_lparen = 0;
static const scg::symbol *S_rparen = 0;
static const scg::symbol *S_eq = 0;
static const scg::symbol *S_leq = 0;
static const scg::symbol *S_geq = 0;
static const scg::symbol *S_neq = 0;
static const scg::symbol *S_lt = 0;
static const scg::symbol *S_gt = 0;

#define GETSYM(s, name) do { if (!(s = scg::symdb_findsymbol_name(name))) { close(sock); return CACLIENT_LIBERROR; } } while (0)

static int sock = -1;

static int errval = CASREP_OK;
int casclient_geterror() { return errval; }
static void seterr(int e) { errval = e; }

const char *
casclient_replystring() {
    int e = casclient_geterror();
	switch (e) {
	case CASREP_OK: return "OK";
	case CASREP_EXPR: return "MathML expression";
	case CASREP_HUP: return "Server hangup";
	case CASREP_MSGERR: return "Bad request";
	case CASREP_SYSERR: return "Server error";
	case CASREP_BADEXPR: return "Malformed expression";
	case CASREP_TOOLARGE: return "Expressions too large";
	case CASREP_INVALIDCMD: return "Command not applicable";
	case CASREP_CONVERR: return "Could not convert Sage output to expression tree";
	case CACLIENT_NETERROR: return "Network error";
	case CASREP_EMPTYMSG: return "Trying to send empty message";
	default: return "(unknown reply code)";
	}
}

int
casclient_connect(void) {
	struct hostent *host;
	struct sockaddr sa;

	int e;

#ifdef WIN32
	WORD vers = MAKEWORD(1,1);
	WSADATA data;
	e = WSAStartup(vers, &data);
	if (e != 0) return CACLIENT_NETERROR;
#endif
	sock = socket(PF_INET, SOCK_STREAM, IPPROTO_TCP);
	if (sock < 0) {
		return CACLIENT_NETERROR;
	}

	host = gethostbyname(CASERVER);
	if (!host) {
		close(sock);
		return CACLIENT_NETERROR;
	}

	sa.sa_family = AF_INET;
	*((unsigned short *)sa.sa_data) = htons(CASERVER_PORT);
	memcpy((void *)(sa.sa_data + sizeof(unsigned short)),
	       host->h_addr_list[0], 4);
	if (connect(sock, &sa, sizeof(sa)) < 0) {
		close(sock);
		return CACLIENT_NETERROR;
	}
	/*if (fcntl(sock, F_SETFL, O_NONBLOCK) < 0) {
		close(sock);
		return CACLIENT_NETERROR;
	}*/

	e = scg::InitializeRecognizer();
	if (FAILURE(e)) {
		close(sock);
		return CACLIENT_LIBERROR;
	}
	GETSYM(S_infin, "infin");
	GETSYM(S_lparen, "lparen");
	GETSYM(S_rparen, "rparen");
	GETSYM(S_plus, "plus");
	GETSYM(S_horzline, "horzline");
	GETSYM(S_eq, "eq");
	GETSYM(S_neq, "neq");
	GETSYM(S_leq, "leq");
	GETSYM(S_geq, "geq");
	GETSYM(S_lt, "lt");
	GETSYM(S_gt, "gt");
	return 0;
}

int
casclient_disconnect(void) {
	if (sock >= 0) {
		int32_t cmd = CASCMD_HUP;
		send(sock, (char *)&cmd, sizeof(cmd), 0);
#ifdef WIN32
		shutdown(sock, 2); // hopefully, SHUT_RDWR == 2
#else
		close(sock); // NB: this doesn't work on win32
#endif
		sock = -1;
	}
	scg::ShutdownRecognizer();
#ifdef WIN32
	WSACleanup();
#endif
	return 0;
}

bool
casclient_isconnected() {
	if (sock < 0) return false;
	struct sockaddr addr;
#ifdef WIN32
	int len;
#else
	socklen_t len;
#endif
	len = sizeof(addr);
	int e = getpeername(sock, &addr, &len);
	if (e < 0) {
		VERBOSE(*verb_out << "caclient: getpeername() returned " << e << " with errno " << errno << " or " << strerror(errno) << std::endl);
	}
	// check error codes
	return (e < 0) ? false : true;
}

std::set<CASOperation> casclient_getoperations(const scg::ExpressionTree *expr)
{
    std::set<CASOperation> ops;
    getOperations(expr, ops);
    return ops;
}

// apply the given operation on the given expression tree
scg::ExpressionTree *
casclient_dooperation(const CASOperation &op, scg::ExpressionTree *expr)
{
    // if we're not connected, try to connect
    if (!casclient_isconnected())
    {
        if (casclient_connect() != 0)
        {
            seterr(CACLIENT_NETERROR);
            return NULL; // cannot connect
        }
    }
    
    const scg::ExpressionTree *varExpr = getVariableExpr(expr, op.variable);

    switch (op.commandCode) {
        case CASCMD_SIMPLIFY:
            return casclient_simplify(expr);
        case CASCMD_EVAL:
            return casclient_evaluate(expr);
        case CASCMD_EVALNUM:
            return casclient_evalnum(expr);
        case CASCMD_EXPAND:
            return casclient_expand(expr);
        case CASCMD_FACTOR:
            return casclient_factor(expr);
        case CASCMD_SOLVE:
            return casclient_solve(expr, varExpr);
        case CASCMD_SOLVESYSTEM:
            return casclient_solvesystem(expr);
        case CASCMD_DETERMINANT:
            return casclient_determinant(expr);
        case CASCMD_INVERSE:
            return casclient_inverse(expr);
        case CASCMD_RANK:
            return casclient_rank(expr);
        case CASCMD_NULLSPACE:
            return casclient_nullspace(expr);
    }
    return NULL;
}

int
casclient_pushexpr(const scg::ExpressionTree *scgtree) {
	int e;
	int32_t rep;
	exprtree *tree = mkexprtree(scgtree);
	if (!tree) {
		return CASREP_SYSERR;
	}
	std::vector<char> buf;
	e = tree->writestream(buf);
	if (e != 0) {
		return e;
	}
	VERBOSE(*verb_out << "caclient: pushing " << buf.size() << "-byte tree " << scgtree->str() << std::endl);
	if (buf.empty()) {
		return CASREP_MSGERR;
	}
	
	rep = CASCMD_EXPR;
	if ((e = fullwrite(sock, (void *)&rep, sizeof(rep))) != sizeof(rep)) {
		return e;
	}
	VERBOSE(*verb_out << "caclient: wrote command code " << rep << std::endl);
	uint32_t wrlen32 = buf.size();
	if ((e = fullwrite(sock, (void *)&wrlen32, sizeof(wrlen32))) != sizeof(wrlen32)) {
		return e;
	}
	VERBOSE(*verb_out << "caclient: wrote length " << buf.size() << std::endl);
	if ((e = fullwrite(sock, (void *)&buf[0], buf.size())) != buf.size()) {
		return e;
	}
	VERBOSE(*verb_out << "caclient: wrote tree buffer" << std::endl);
	if ((e = fullread(sock, (void *)&rep, sizeof(rep))) != sizeof(rep)) {
		return e;
	}
	VERBOSE(*verb_out << "caclient: got reply " << rep << std::endl);
	return rep;
}

scg::ExpressionTree *
exprcmd(int32_t cmd) {
	VERBOSE(*verb_out << "caclient: executing command " << cmd << std::endl);
	int e;
	int32_t rep;
	if ((e = fullwrite(sock, (void *)&cmd, sizeof(cmd))) != sizeof(cmd)) {
		seterr(e);
		return 0;
	}
	VERBOSE(*verb_out << "caclient: wrote command code\n");
	if ((e = fullread(sock, (void *)&rep, sizeof(rep))) != sizeof(rep)) {
		seterr(e);
		return 0;
	}
	VERBOSE(*verb_out << "caclient: read response " << rep << std::endl);
	if (rep == CASREP_EXPR) {
		exprtree *tree;
		if ((e = readexpr(sock, &tree)) != CASREP_OK) {
			seterr(e);
			return 0;
		}
		scg::ExpressionTree *expr = tree->mkscgtree();
		if (!expr) {
			seterr(CASREP_SYSERR);
		}
		else {
			VERBOSE(*verb_out << "caclient: read tree " << expr->str() << std::endl);
			seterr(CASREP_OK);
		}
		return expr;
	}
	else {
		seterr(rep);
		return 0;
	}
}

scg::ExpressionTree *
casclient_evaluate(scg::ExpressionTree *expr) {
	seterr(casclient_pushexpr(expr));
	if (casclient_geterror() != CASREP_OK) return 0;
	return exprcmd(CASCMD_EVAL);
}

scg::ExpressionTree *
casclient_evalnum(scg::ExpressionTree *expr) {
	seterr(casclient_pushexpr(expr));
	if (casclient_geterror() != CASREP_OK) return 0;
	return exprcmd(CASCMD_EVALNUM);
}

scg::ExpressionTree *
casclient_simplify(scg::ExpressionTree *expr) {
	seterr(casclient_pushexpr(expr));
	if (casclient_geterror() != CASREP_OK) return 0;
	return exprcmd(CASCMD_SIMPLIFY);
}

scg::ExpressionTree *
casclient_expand(scg::ExpressionTree *expr) {
	seterr(casclient_pushexpr(expr));
	if (casclient_geterror() != CASREP_OK) return 0;
	return exprcmd(CASCMD_EXPAND);
}

scg::ExpressionTree *
casclient_factor(scg::ExpressionTree *expr) {
	seterr(casclient_pushexpr(expr));
	if (casclient_geterror() != CASREP_OK) return 0;
	return exprcmd(CASCMD_FACTOR);
}

scg::ExpressionTree *
casclient_solve(scg::ExpressionTree *expr, const scg::ExpressionTree *var) {
	seterr(casclient_pushexpr(expr));
	if (casclient_geterror() != CASREP_OK) return 0;
	seterr(casclient_pushexpr(var));
	if (casclient_geterror() != CASREP_OK) return 0;
	return exprcmd(CASCMD_SOLVE);
}

scg::ExpressionTree *
casclient_solvesystem(scg::ExpressionTree *expr) {
 	seterr(casclient_pushexpr(expr));
 	if (casclient_geterror() != CASREP_OK) return 0;
 	return exprcmd(CASCMD_SOLVESYSTEM);
}

scg::ExpressionTree *
casclient_determinant(scg::ExpressionTree *expr) {
	seterr(casclient_pushexpr(expr));
	if (casclient_geterror() != CASREP_OK) return 0;
	return exprcmd(CASCMD_DETERMINANT);
}

scg::ExpressionTree *
casclient_inverse(scg::ExpressionTree *expr) {
	seterr(casclient_pushexpr(expr));
	if (casclient_geterror() != CASREP_OK) return 0;
	return exprcmd(CASCMD_INVERSE);
}

scg::ExpressionTree *
casclient_rank(scg::ExpressionTree *expr) {
	seterr(casclient_pushexpr(expr));
	if (casclient_geterror() != CASREP_OK) return 0;
	return exprcmd(CASCMD_RANK);
}

scg::ExpressionTree *
casclient_nullspace(scg::ExpressionTree *expr) {
	seterr(casclient_pushexpr(expr));
	if (casclient_geterror() != CASREP_OK) return 0;
	return exprcmd(CASCMD_NULLSPACE);
}
