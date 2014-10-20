#include "profman.h"

#include <algorithm>
#include <fstream>
#include <functional>
#include <sstream>
#include <string>
#include <iostream>

#include "confusion.h"
#include "debug.h"
#include "feat.h"
#include "memory.h"
#include "profcache.h"
#include "recog.h"
#include "symtab.h"
#include "vector-io.h"
#include "utils.h"


namespace scg
{


ProfileManager::ProfileManager(const SymbolTable &symtab)
{
    for (SymbolTable::const_iterator i = symtab.begin(); i != symtab.end(); ++i) {
        symbol_table.insert(std::make_pair(i->unicode_char(), *i));
    }
}


//ProfileManager::~ProfileManager()
//{
    //for (std::map<PrototypeId, PrototypeData>::iterator i = prototypes.begin(); i != prototypes.end(); ++i) {
    //    delete i->second.correlations;
    //}
//}


int
ProfileManager::ImportPrototypeCorrelations(const std::vector<InterProfileCorrelation> &correlations)
{
   for (std::vector<InterProfileCorrelation>::const_iterator i = correlations.begin(); i != correlations.end(); ++i) {
      //prototypes[i->from_id].correlations->push_back(PrototypeCorrelation(i->to_id, i->correlation));
		Prototype *proto = GetPrototype(i->from_id);
		if (proto) {
        	proto->correlations.push_back(PrototypeCorrelation(i->to_id, i->correlation));
		}
   }
   return 0;
}


/*
int
ProfileManager::ImportPrototypeConfusion(const std::vector<InterProfileCorrelation> *confusions, unsigned n)
{
    for (unsigned j = 0; j < n; j++) {
        for (std::vector<InterProfileCorrelation>::const_iterator i = confusions[j].begin(); i != confusions[j].end(); ++i) {
            GetPrototype(i->from_id)->confusions[j].push_back(PrototypeCorrelation(i->to_id, i->correlation));
        }
    }
    return 0;
}
*/

int
ProfileManager::UpdateProfileCorrelations(const Profile &from, const Profile &to)
{
    int e = E_INVALID;
    std::vector<InterProfileCorrelation> correlations;
    if (AreProfileCorrelationsCached(from.Id(), to.Id())) {
        e = ImportProfileCorrelations(from.Id(), to.Id(), correlations);
    }
    if (FAILURE(e)) {
        e = ComputeProfileCorrelations(from, to, correlations);
        if (FAILURE(e)) {
            return e;
        }
        e = CacheProfileCorrelations(*this, from.Id(), to.Id(), correlations);
        if (FAILURE(e)) {
            return e;
        }
    }
    
    return ImportPrototypeCorrelations(correlations);
}


/*
int
ProfileManager::UpdateProfileConfusion(const Profile &from, const Profile &to)
{
    int e = E_INVALID;
    std::vector<InterProfileCorrelation> *confusions = DEBUG_NEW std::vector<InterProfileCorrelation>[NumMatchers];
    if (IsProfileConfusionCached(from.Id(), to.Id())) {
        e = ImportProfileConfusion(from.Id(), to.Id(), confusions, NumMatchers);
    }
    if (FAILURE(e)) {
        e = ComputeProfileConfusion(*this, from, to, confusions);
        if (FAILURE(e)) {
            return e;
        }
        e = CacheProfileConfusion(from.Id(), to.Id(), confusions, NumMatchers);
        if (FAILURE(e)) {
            return e;
        }
    }
    
    e = ImportPrototypeConfusion(confusions, NumMatchers);
    
    delete[] confusions;
    
    return e;
}
*/

int
ProfileManager::AddProfile(Profile &prof)
{
    profiles.push_back(prof);
    
    const Profile &profile = profiles.back();
    
    // Update prototype map
    /*for (Profile::const_iterator symbol = profile.begin(); symbol != profile.end(); ++symbol) {
        for (ProfileSymbol::const_iterator prototype = symbol->begin(); prototype != symbol->end(); ++prototype) {
            PrototypeData pd;
            pd.prototype = &(*prototype);
            pd.correlations = DEBUG_NEW std::vector<PrototypeCorrelation>;
            DEBUG_ONLY(debug_out << "importing " << prototype->id << std::endl);
                        
            std::pair<std::map<PrototypeId, PrototypeData>::iterator, bool> p = prototypes.insert(std::make_pair(prototype->id, pd));
            if (p.second == false) {
                return E_INTERNAL;
            }
        }
    }*/
/*
    for (Profile::const_iterator symbol1 = profile.begin(); symbol1 != profile.end(); ++symbol1) {
        for (ProfileSymbol::const_iterator prototype1 = symbol1->begin(); prototype1 != symbol1->end(); ++prototype1) {
            std::map<PrototypeId, PrototypeData>::const_iterator i = prototypes.find(prototype1->id);
            if (i == prototypes.end()) {
                DEBUG_ONLY(debug_out << "lost " << prototype1->id << std::endl);
            }
            else {
                const PrototypeData &pd1 = i->second;
                assert(pd1.prototype == &(*prototype1));
            }
        }
    }
*/
        
    // Update feature correlation
    //for (ProfileContainer::const_iterator p = profiles.begin(); p != profiles.end(); ++p) {
	 for (unsigned i = 0; i < profiles.size(); i++) {
        int e = UpdateProfileCorrelations(profiles[i], profile);
        if (FAILURE(e)) {
            return e;
        }
        
        if (profiles[i] != profile) {
            e = UpdateProfileCorrelations(profile, profiles[i]);
            if (FAILURE(e)) {
                return e;
            }
        }
    }

    // Update confusion symbols
    /*for (ProfileContainer::const_iterator p = profiles.begin(); p != profiles.end(); ++p) {
        int e = UpdateProfileConfusion(*p, profile);
        if (FAILURE(e)) {
            return e;
        }
        
        if (*p != profile) {
            e = UpdateProfileConfusion(profile, *p);
            if (FAILURE(e)) {
                return e;
            }
        }
    }*/
    
    return 0;
}


SymbolClass
ProfileManager::GetSymbolClassByUnicodeChar(wchar_t unicode) const
{
    return SymbolClass(SymbolClass::const_iterator(begin(), end(), unicode),
                       SymbolClass::const_iterator(begin(), end(), unicode));
}


Prototype *
ProfileManager::GetPrototype(const PrototypeId &id)
{
    for (ProfileContainer::iterator i = profiles.begin(); i != profiles.end(); ++i) {
        if (i->Id() == id.profile) {
            ProfileSymbol &symbol = *(i->begin() + id.entry);
            return &(*(symbol.begin() + id.prototype));
        }
    }

    return 0;    
}

const Prototype *
ProfileManager::GetPrototypeById(const PrototypeId &id) const
{
    //std::map<PrototypeId, PrototypeData>::const_iterator i;
    //i = prototypes.find(id);
    //if (i == prototypes.end()) {
    //    return 0;
    //}
    
    //return i->second.prototype;

    for (ProfileContainer::const_iterator i = profiles.begin(); i != profiles.end(); ++i) {
        if (i->Id() == id.profile) {
            const ProfileSymbol &symbol = *(i->begin() + id.entry);
            return &(*(symbol.begin() + id.prototype));
        }
    }

    return 0;
}


const ProfileSymbol *
ProfileManager::GetSymbolByName(const std::string &name) const
{
    for (const_iterator i = begin(); i != end(); ++i) {
        if (GetInfoBySymbol(*i).name() == name) {
            return &(*i);
        }
    }
    return 0;
}

const SymbolInfo &
ProfileManager::GetInfoByPrototypeId(const PrototypeId &id) const
{
    const Prototype *prototype = GetPrototypeById(id);
    if (!prototype) {
        //DEBUG_ONLY(debug_out << "missing " << id << std::endl);
        throw E_NOTFOUND;
    }
    return GetInfoByPrototype(*prototype);
}


const SymbolInfo &
ProfileManager::GetInfoByPrototype(const Prototype &prototype) const
{
    return GetInfoByUnicodeChar(prototype.symbol->unicode_char());
}

const SymbolInfo &
ProfileManager::GetInfoBySymbol(const ProfileSymbol &symbol) const
{
    return GetInfoByUnicodeChar(symbol.unicode_char());
}

const SymbolInfo &
ProfileManager::GetInfoByUnicodeChar(wchar_t ch) const
{
    std::map<wchar_t, SymbolInfo>::const_iterator i = symbol_table.find(ch);
    if (i == symbol_table.end()) {
        throw E_NOTFOUND;
    }
    else {
        return i->second;
    }
}


PrototypeCorrelations
ProfileManager::GetPrototypeCorrelationsFor(const PrototypeId &id) const
{
    static std::vector<PrototypeCorrelation> EmptyCorrelations;
/*
    std::map<PrototypeId, PrototypeData>::const_iterator i;
    i = prototypes.find(id);
    if (i == prototypes.end()) {
        return PrototypeCorrelations(EmptyCorrelations);
    }
    return PrototypeCorrelations(*i->second.correlations);
    */
    
    const Prototype *prototype = GetPrototypeById(id);
    if (prototype) {
        return PrototypeCorrelations(prototype->correlations);
    }
    return PrototypeCorrelations(EmptyCorrelations);
}

/*
PrototypeCorrelations
ProfileManager::GetPrototypeConfusionsFor(const PrototypeId &id, int matcher) const
{
    /*static std::vector<const Prototype *> EmptyConfusion;
    
    std::map<PrototypeId, PrototypeData>::const_iterator i;
    i = prototypes.find(id);
    if (i == prototypes.end()) {
        return ConfusionPrototypes(EmptyConfusion);
    }
    return ConfusionPrototypes(i->second.confusion_symbols);* /
    static std::vector<PrototypeCorrelation> EmptyConfusions;
    const Prototype *prototype = GetPrototypeById(id);
    if (prototype) {
        return PrototypeCorrelations(prototype->confusions[matcher]);
    }
    return PrototypeCorrelations(EmptyConfusions);
}*/


ProfileManager::iterator
ProfileManager::begin()
  { return iterator(profiles); }


ProfileManager::iterator
ProfileManager::end()
{
    return iterator(profiles, profiles.end());
}

ProfileManager::const_iterator
ProfileManager::begin() const
{
    return const_iterator(profiles);
}


ProfileManager::const_iterator
ProfileManager::end() const
{
    if (profiles.empty()) {
        return begin();
    }
    return const_iterator(profiles, profiles.end());
}

size_t
ProfileManager::size() const
{
    size_t count = 0;
    for (ProfileContainer::const_iterator i = profiles.begin(); i != profiles.end(); ++i) {
        count += i->size();
    }
    return count;
}


bool
ProfileManager::empty() const
{
    if (profiles.empty()) {
        return true;
    }

    return scg::find_if(profiles.begin(), profiles.end(), std::not1(std::mem_fun_ref(&Profile::empty))) == profiles.end();
}



int
CreateDefaultProfileManager(ProfileManager *& profman)
{
    int error = 0;
    
    profman = 0;
    
    error = scg::recognizer_initialize();
    if (FAILURE(error)) {
        return error;
    }
    
    std::string training_path;
    scg::GetTrainingPath(training_path);
    
    scg::SymbolTable symbol_table;
    error = scg::CreateDefaultSymbolTable(symbol_table);
    if (FAILURE(error)) {
        return error;
    }
    
    std::ifstream list_file((training_path + "/profiles2.list").c_str());
    if (!list_file.is_open()) {
        return E_NOTFOUND;
    }
    
    profman = DEBUG_NEW scg::ProfileManager(symbol_table);
    if (!profman) {
        return E_OUTOFMEM;
    }

    static const size_t BufSz = 1024;
    char buf[BufSz];
    
    unsigned n = 0;
    
    while (list_file) {
        std::stringstream filename;
        
        do {
            memset(buf, 0, BufSz);
            list_file.clear(std::ios_base::goodbit);
            list_file.getline(buf, BufSz);
            if (list_file.bad() || list_file.eof()) {
                break;
            }
            filename << buf;
        } while (list_file.fail());
        
        if (list_file.bad()) {
            delete profman;
            return E_IO;
        }
        
        if (!filename.str().empty()) {
            std::ifstream profile_file;
        
            error = OpenProfileFile(filename.str(), profile_file);
            if (FAILURE(error)) {
                DEBUG_ONLY(scg::debug_out << "could not open " << filename.str() << std::endl);
                std::cerr << "could not open " << filename.str() << std::endl;
                continue;
            }
                      
            try {
                scg::Profile profile;
                profile_file >> profile;
					 if (profile.Id() != scg::Profile::BadId) {
	                error = profman->AddProfile(profile);
   	             if (FAILURE(error)) {
      	              throw error;
        	   	    }
         	       n++;
                }
                // For the first non-default profile loaded, load auxiliary data
                // XXX: This is only correct because MathBrush only lets the user
                // select a single profile rather than multiple ones.  Some more
                // general solution is needed.
                /*if (n == 2) {
                    std::stringstream auxdata_filename;
                    std::string::size_type i = filename.find_last_of("/\\");
                    if (i != std::string::npos) {
                        auxdata_filename << filename.substr(0, i);
                    }
                    
                    auxdata_filename << profile.Id() << ".aux";
                    
                    std::ifstream auxdata(auxdata_filename.str().c_str());
                    if (auxdata.is_open()) {
                        for (;;) {
                            std::string dummy_id;
                            auxdata >> dummy_id;
                            if (auxdata.eof()) {
                                break;
                            }
                            
                            std::stringstream idss;
                            idss << auxdata;
                            
                            PrototypeId id;
                            idss >> id;
                            
                            
                        }
                    }
                }*/
            }
            catch (int e) {
                DEBUG_ONLY(scg::debug_out << "failure loading " << filename.str() << "(" << e << ")" << std::endl);
            }
        }
    }
    
    return 0;
}


}

