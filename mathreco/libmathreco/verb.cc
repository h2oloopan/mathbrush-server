#ifndef VERB_H_
#define VERB_H_


#include "dlldecl.h"
#include "error.h"
#include "verb.h"

#include <ostream>
#include <cstdlib>
#include <fstream>


namespace scg
{


std::ostream *verb_out;
static bool delete_verb_out = false;
int verb_level = 0;
static bool verb_init = false;


DLLDECL void
SetVerbosity(int level)
{
    verb_level = level;
}


DLLDECL void
NoVerboseOutput()
{
    if (delete_verb_out) {
        delete verb_out;
        delete_verb_out = false;
    }
    verb_out = 0;
}


DLLDECL void
VerboseOutputToStream(std::ostream &os)
{
    NoVerboseOutput();
    verb_out = &os;
    delete_verb_out = false;
}


DLLDECL void
VerboseOutputToFile(const char *f)
{
    NoVerboseOutput();
    std::ofstream *ofs = new std::ofstream(f);
    if (!ofs->is_open() || ofs->bad()) {
        delete ofs;
        throw E_IO;
    }
	verb_out = ofs;
    delete_verb_out = true;
    if (!verb_init) {
        std::atexit(&NoVerboseOutput);
        verb_init = true;
    }
}



}

#endif
