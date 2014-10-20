#include <fstream>
#include <string>

#include "error.h"
#include "info.h"
#include "profile.h"
#include "symtab.h"
#include "vector-io.h"
#include "utils.h"


namespace scg
{


int
CreateDefaultSymbolTable(SymbolTable &symtab)
{
    std::string training_path;
    
    int err = scg::GetTrainingPath(training_path);
    if (FAILURE(err)) {
        return err;
    }
    
    std::ifstream desc_file((training_path + "/basic-symbols.desc").c_str());
    try {
        desc_file >> symtab;
    }
    catch (int e) {
        return e;
    }
    
    return 0;
}


}

