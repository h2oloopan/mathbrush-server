#include "parms.h"

#include <cstdlib>
#include <fstream>
#include <map>
#include <vector>
#include <string>

#include "memory.h"
#include "error.h"
#include "utils.h"


namespace scg
{

static std::map<std::string, int> *int_parameters_min = 0;
static std::map<std::string, unsigned> *unsigned_parameters_min = 0;
static std::map<std::string, double> *double_parameters_min = 0;

static std::map<std::string, int> *int_parameters_max = 0;
static std::map<std::string, unsigned> *unsigned_parameters_max = 0;
static std::map<std::string, double> *double_parameters_max = 0;


static std::map<std::string, int> *int_parameters = 0;
static std::map<std::string, unsigned> *unsigned_parameters = 0;
static std::map<std::string, double> *double_parameters = 0;
static std::map<std::string, std::string> *string_parameters = 0;
static bool parameters_loaded = false;


static int
SetParameterMin(const std::string &name, int value)
{
    if (!int_parameters_min) {
        int_parameters_min = new std::map<std::string, int>;
    }
    
    (*int_parameters_min)[name] = value;
        
    return 0;
}

static int
SetParameterMin(const std::string &name, unsigned value)
{
    if (!unsigned_parameters_min) {
        unsigned_parameters_min = new std::map<std::string, unsigned>;
    }
    
    (*unsigned_parameters_min)[name] = value;
        
    return 0;
}

static int
SetParameterMin(const std::string &name, double value)
{
    if (!double_parameters_min) {
        double_parameters_min = new std::map<std::string, double>;
    }

    (*double_parameters_min)[name] = value;
        
    return 0;
}


static int
SetParameterMax(const std::string &name, int value)
{
    if (!int_parameters_max) {
        int_parameters_max = new std::map<std::string, int>;
    }
    
    (*int_parameters_max)[name] = value;
        
    return 0;
}

static int
SetParameterMax(const std::string &name, unsigned value)
{
    if (!unsigned_parameters_max) {
        unsigned_parameters_max = new std::map<std::string, unsigned>;
    }
    
    (*unsigned_parameters_max)[name] = value;
        
    return 0;
}

static int
SetParameterMax(const std::string &name, double value)
{
    if (!double_parameters_max) {
        double_parameters_max = new std::map<std::string, double>;
    }

    (*double_parameters_max)[name] = value;
        
    return 0;
}


static int
GetParmPath(std::string &path)
{
    int error;
    error = GetProfilePath(path);
    if (FAILURE(error)) {
        return error;
    }
    
    path += "/parms2.list";
    
    return 0;
}

template <typename T>
static void
reset(T *& p) {
	delete p;
	p = 0;
}

int
ReleaseParameters()
{
    reset(int_parameters);
    reset(int_parameters_min);
    reset(int_parameters_max);
    reset(unsigned_parameters);
    reset(unsigned_parameters_min);
    reset(unsigned_parameters_max);
    reset(double_parameters);
    reset(double_parameters_min);
    reset(double_parameters_max);
    reset(string_parameters);
    return 0;
}

int
RestoreParameters()
{
    std::string parm_path;
    int error;
    
    error = GetParmPath(parm_path);
    if (FAILURE(error)) {
        return error;
    }

    std::ifstream ifs(parm_path.c_str());
    if (!ifs.is_open()) {
        return E_NOTFOUND;
    }
    
    for (;;) {
        std::string type;
        ifs >> type;
        if (ifs.eof()) {
            break;
        }
        
        if (type == "//") {
            std::string dummy;
            std::getline(ifs, dummy);
            continue;
        }
        
        std::string name;
        ifs >> name;
        
        if (type == "int") {
            int val, minval, maxval;
            ifs >> val >> minval >> maxval;
            if (ifs.eof() || ifs.bad()) {
                return E_IO;
            }
            SetParameter(name, val, false);
            SetParameterMin(name, minval);
            SetParameterMax(name, maxval);
        }
        else if (type == "unsigned") {
            unsigned val, minval, maxval;
            ifs >> val >> minval >> maxval;
            if (ifs.eof() || ifs.bad()) {
                return E_IO;
            }
            SetParameter(name, val, false);
            SetParameterMin(name, minval);
            SetParameterMax(name, maxval);
        }
        else if (type == "double") {
            double val, minval, maxval;
            ifs >> val >> minval >> maxval;
            if (ifs.eof() || ifs.bad()) {
                return E_IO;
            }
            SetParameter(name, val, false);
            SetParameterMin(name, minval);
            SetParameterMax(name, maxval);
        }
        else if (type == "string") {
            std::string val;
            ifs >> val;
            if (ifs.eof() || ifs.bad()) {
                return E_IO;
            }
            SetParameter(name, val, false);        
        }
    }
    
    parameters_loaded = true;
    
//    _onexit(&ReleaseParameters);
    
    return 0;
}


static int
PersistParameters()
{
    std::string parm_path;
    int error;
    
    error = GetParmPath(parm_path);
    if (FAILURE(error)) {
        return error;
    }
    
    std::ofstream ofs(parm_path.c_str());
    if (!ofs.is_open()) {
        return E_IO;
    }
    
    if (int_parameters) {
        for (std::map<std::string, int>::const_iterator i = int_parameters->begin(); i != int_parameters->end(); ++i) {
            ofs << "int " << i->first << " " << i->second << std::endl;
        }
    }
    if (unsigned_parameters) {
        for (std::map<std::string, unsigned>::const_iterator i = unsigned_parameters->begin(); i != unsigned_parameters->end(); ++i) {
            ofs << "unsigned " << i->first << " " << i->second << std::endl;
        }
    }
    if (double_parameters) {
        for (std::map<std::string, double>::const_iterator i = double_parameters->begin(); i != double_parameters->end(); ++i) {
            ofs << "double " << i->first << " " << i->second << std::endl;
        }
    }
    if (string_parameters) {
        for (std::map<std::string, std::string>::const_iterator i = string_parameters->begin(); i != string_parameters->end(); ++i) {
            ofs << "string " << i->first << " " << i->second << std::endl;
        }
    }
        
    return 0;
}


int
GetParameterInt(const std::string &name)
{
    if (!parameters_loaded) {
        RestoreParameters();
    }
    
    if (!int_parameters) {
        throw E_NOTFOUND;
    }
    
    std::map<std::string, int>::iterator parm = int_parameters->find(name);
    if (parm == int_parameters->end()) {
        throw E_NOTFOUND;
    }
    return parm->second;
}

unsigned
GetParameterUnsigned(const std::string &name)
{
    if (!parameters_loaded) {
        RestoreParameters();
    }
    
    if (!unsigned_parameters) {
        throw E_NOTFOUND;
    }
    
    std::map<std::string, unsigned>::iterator parm = unsigned_parameters->find(name);
    if (parm == unsigned_parameters->end()) {
        throw E_NOTFOUND;
    }
    return parm->second;
}

double
GetParameterDouble(const std::string &name)
{
    if (!parameters_loaded) {
        RestoreParameters();
    }
    
    if (!double_parameters) {
        throw E_NOTFOUND;
    }
    
    std::map<std::string, double>::iterator parm = double_parameters->find(name);
    if (parm == double_parameters->end()) {
        throw E_NOTFOUND;
    }
    return parm->second;
}

std::string
GetParameterString(const std::string &name)
{
    if (!parameters_loaded) {
        RestoreParameters();
    }

    if (!string_parameters) {
        throw E_NOTFOUND;
    }
        
    std::map<std::string, std::string>::iterator parm = string_parameters->find(name);
    if (parm == string_parameters->end()) {
        throw E_NOTFOUND;
    }
    return parm->second;
}


int
GetParameterIntMin(const std::string &name)
{
    if (!parameters_loaded) {
        RestoreParameters();
    }
    
    if (!int_parameters_min) {
        throw E_NOTFOUND;
    }
    
    std::map<std::string, int>::iterator parm = int_parameters_min->find(name);
    if (parm == int_parameters_min->end()) {
        throw E_NOTFOUND;
    }
    return parm->second;
}

unsigned
GetParameterUnsignedMin(const std::string &name)
{
    if (!parameters_loaded) {
        RestoreParameters();
    }
    
    if (!unsigned_parameters_min) {
        throw E_NOTFOUND;
    }
    
    std::map<std::string, unsigned>::iterator parm = unsigned_parameters_min->find(name);
    if (parm == unsigned_parameters_min->end()) {
        throw E_NOTFOUND;
    }
    return parm->second;
}

double
GetParameterDoubleMin(const std::string &name)
{
    if (!parameters_loaded) {
        RestoreParameters();
    }
    
    if (!double_parameters_min) {
        throw E_NOTFOUND;
    }
    
    std::map<std::string, double>::iterator parm = double_parameters_min->find(name);
    if (parm == double_parameters_min->end()) {
        throw E_NOTFOUND;
    }
    return parm->second;
}


int
GetParameterIntMax(const std::string &name)
{
    if (!parameters_loaded) {
        RestoreParameters();
    }
    
    if (!int_parameters_max) {
        throw E_NOTFOUND;
    }
    
    std::map<std::string, int>::iterator parm = int_parameters_max->find(name);
    if (parm == int_parameters_max->end()) {
        throw E_NOTFOUND;
    }
    return parm->second;
}

unsigned
GetParameterUnsignedMax(const std::string &name)
{
    if (!parameters_loaded) {
        RestoreParameters();
    }
    
    if (!unsigned_parameters_max) {
        throw E_NOTFOUND;
    }
    
    std::map<std::string, unsigned>::iterator parm = unsigned_parameters_max->find(name);
    if (parm == unsigned_parameters_max->end()) {
        throw E_NOTFOUND;
    }
    return parm->second;
}

double
GetParameterDoubleMax(const std::string &name)
{
    if (!parameters_loaded) {
        RestoreParameters();
    }
    
    if (!double_parameters_max) {
        throw E_NOTFOUND;
    }
    
    std::map<std::string, double>::iterator parm = double_parameters_max->find(name);
    if (parm == double_parameters_max->end()) {
        throw E_NOTFOUND;
    }
    return parm->second;
}


int
SetParameter(const std::string &name, int value, bool persist)
{
    if (!int_parameters) {
        int_parameters = new std::map<std::string, int>;
    }
    
    (*int_parameters)[name] = value;
        
    if (persist) {
        PersistParameters();
    }
    
    return 0;
}

int
SetParameter(const std::string &name, unsigned value, bool persist)
{
    if (!unsigned_parameters) {
        unsigned_parameters = new std::map<std::string, unsigned>;
    }
    
    (*unsigned_parameters)[name] = value;
        
    if (persist) {
        PersistParameters();
    }
    
    return 0;
}

int
SetParameter(const std::string &name, double value, bool persist)
{
    if (!double_parameters) {
        double_parameters = new std::map<std::string, double>;
    }

    (*double_parameters)[name] = value;
        
    if (persist) {
        PersistParameters();
    }
    
    return 0;
}

int
SetParameter(const std::string &name, const std::string &value, bool persist)
{
    if (!string_parameters) {
        string_parameters = new std::map<std::string, std::string>;
    }

    (*string_parameters)[name] = value;
        
    if (persist) {
        PersistParameters();
    }
    
    return 0;
}


const std::map<std::string, int> &
GetIntParameters()
{
    if (!parameters_loaded) {
        RestoreParameters();
    }
    
    if (!int_parameters) {
        int_parameters = new std::map<std::string, int>;        
    }
    
    return *int_parameters;
}

const std::map<std::string, unsigned> &
GetUnsignedParameters()
{
    if (!parameters_loaded) {
        RestoreParameters();
    }
    
    if (!unsigned_parameters) {
        unsigned_parameters = new std::map<std::string, unsigned>;        
    }
    
    return *unsigned_parameters;
}

const std::map<std::string, double> &
GetDoubleParameters()
{
    if (!parameters_loaded) {
        RestoreParameters();
    }
    
    if (!double_parameters) {
        double_parameters = new std::map<std::string, double>;        
    }
    
    return *double_parameters;
}

const std::map<std::string, std::string> &
GetStringParameters()
{
    if (!parameters_loaded) {
        RestoreParameters();
    }
    
    if (!string_parameters) {
        string_parameters = new std::map<std::string, std::string>;        
    }
    
    return *string_parameters;
}


const double UNINITIALIZED_PARM_DOUBLE = 0.0;
const int UNINITIALIZED_PARM_INT = 0;
const unsigned UNINITIALIZED_PARM_UNSIGNED = 0;
const std::string UNINITIALIZED_PARM_STRING = "";

static std::vector<std::pair<std::string, int *> > &
int_parms() {
	static std::vector<std::pair<std::string, int *> > int_parms_;
	return int_parms_;
}
static std::vector<std::pair<std::string, unsigned *> > &
unsigned_parms() {
	static std::vector<std::pair<std::string, unsigned *> > unsigned_parms_;
	return unsigned_parms_;
}
static std::vector<std::pair<std::string, double *> > &
double_parms() {
	static std::vector<std::pair<std::string, double *> > double_parms_;
	return double_parms_;
}
static std::vector<std::pair<std::string, std::string *> > &
string_parms() {
	static std::vector<std::pair<std::string, std::string *> > string_parms_;
	return string_parms_;
}

int
RegisterParameterInt(const std::string &name, int *p)
{
	int_parms().push_back(std::make_pair(name, p));
	return UNINITIALIZED_PARM_INT;
}

unsigned
RegisterParameterUnsigned(const std::string &name, unsigned *p)
{
	unsigned_parms().push_back(std::make_pair(name, p));
	return UNINITIALIZED_PARM_UNSIGNED;
}

double
RegisterParameterDouble(const std::string &name, double *p)
{
	double_parms().push_back(std::make_pair(name, p));
	return UNINITIALIZED_PARM_DOUBLE;
}

std::string
RegisterParameterString(const std::string &name, std::string *p)
{
	string_parms().push_back(std::make_pair(name, p));
	return UNINITIALIZED_PARM_STRING;
}


std::vector<parm_cb> &
callbacks() {
	static std::vector<parm_cb> callbacks_;
	return callbacks_;
}

int
RegisterParameterCallback(parm_cb f)
{
	callbacks().push_back(f);
	return 0;
}


int
RestoreRegisteredParms() {
	int e = RestoreParameters();
	if (FAILURE(e)) return e;
	for (std::vector<std::pair<std::string, int *> >::iterator i = int_parms().begin(); i != int_parms().end(); ++i) {
		*(i->second) = GetParameterInt(i->first);
	}
	for (std::vector<std::pair<std::string, unsigned *> >::iterator i = unsigned_parms().begin(); i != unsigned_parms().end(); ++i) {
		*(i->second) = GetParameterUnsigned(i->first);
	}
	for (std::vector<std::pair<std::string, double *> >::iterator i = double_parms().begin(); i != double_parms().end(); ++i) {
		*(i->second) = GetParameterDouble(i->first);
	}
	for (std::vector<std::pair<std::string, std::string *> >::iterator i = string_parms().begin(); i != string_parms().end(); ++i) {
		*(i->second) = GetParameterString(i->first);
	}
	for (std::vector<parm_cb>::iterator i = callbacks().begin(); i != callbacks().end(); ++i) {
		(*i)();
	}
	return 0;
}



}

