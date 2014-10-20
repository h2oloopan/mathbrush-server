#ifndef PARMS_H_
#define PARMS_H_


#include <map>
#include <string>


namespace scg
{


const std::map<std::string, int> &GetIntParameters();
const std::map<std::string, unsigned> &GetUnsignedParameters();
const std::map<std::string, double> &GetDoubleParameters();
const std::map<std::string, std::string> &GetStringParameters();

int GetParameterInt(const std::string &name);
unsigned GetParameterUnsigned(const std::string &name);
double GetParameterDouble(const std::string &name);
std::string GetParameterString(const std::string &name);

int GetParameterIntMin(const std::string &name);
unsigned GetParameterUnsignedMin(const std::string &name);
double GetParameterDoubleMin(const std::string &name);
std::string GetParameterStringMin(const std::string &name);

int GetParameterIntMax(const std::string &name);
unsigned GetParameterUnsignedMax(const std::string &name);
double GetParameterDoubleMax(const std::string &name);
std::string GetParameterStringMax(const std::string &name);


int SetParameter(const std::string &name, int value, bool persist = true);
int SetParameter(const std::string &name, unsigned value, bool persist = true);
int SetParameter(const std::string &name, double value, bool persist = true);
int SetParameter(const std::string &name, const std::string &value, bool persist = true);


extern const double UNINITIALIZED_PARM_DOUBLE;
extern const int UNINITIALIZED_PARM_INT;
extern const unsigned UNINITIALIZED_PARM_UNSIGNED;
extern const std::string UNINITIALIZED_PARM_STRING;

int RegisterParameterInt(const std::string &name, int *p);
unsigned RegisterParameterUnsigned(const std::string &name, unsigned *p);
double RegisterParameterDouble(const std::string &name, double *p);
std::string RegisterParameterDouble(const std::string &name, std::string *p);

int RestoreRegisteredParms();
int RestoreParameters();
int ReleaseParameters();


typedef void (*parm_cb)();
int RegisterParameterCallback(parm_cb f);

}


#endif

