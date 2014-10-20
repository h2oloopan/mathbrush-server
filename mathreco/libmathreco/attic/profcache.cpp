#include "profcache.h"

#include <fstream>
#include <sstream>

#include "confusion.h"
#include "error.h"
#include "feat.h"
#include "vector-io.h"
#include "utils.h"
#include "md5.h"


namespace scg
{


static std::istream &
operator>>(std::istream &is, scg::InterProfileCorrelation &correlation)
{
    is >> correlation.from_id >> correlation.to_id >> correlation.correlation;
    if (!is) {
        throw E_IO;
    }
    
    return is;
}

static std::ostream &
operator<<(std::ostream &os, const scg::InterProfileCorrelation &correlation)
{
    os << correlation.from_id << " " << correlation.to_id << " " << correlation.correlation;
    if (!os) {
        throw E_IO;
    }
    
    return os;
}


static std::istream &
operator>>(std::istream &is, scg::InterProfileConfusion &confusion)
{
    is >> confusion.from_id >> confusion.to_id;
    if (!is) {
        throw E_IO;
    }
    
    return is;
}

static std::ostream &
operator<<(std::ostream &os, const scg::InterProfileConfusion &confusion)
{
    os << confusion.from_id << " " << confusion.to_id;
    if (!os) {
        throw E_IO;
    }
    
    return os;
}


static bool cache_init = false;

typedef std::map<std::pair<ProfileId, ProfileId>, std::string> CorrelationCache;
static CorrelationCache correlation_cache;

typedef std::map<std::pair<ProfileId, ProfileId>, std::string> ConfusionCache;
static ConfusionCache confusion_cache;


static bool
CacheInitialized()
{
    return cache_init;
}

static int
WriteCacheSpecFile()
{
    int e;
    std::string training_path;
    e = GetTrainingPath(training_path);
    if (FAILURE(e)) {
        return e;
    }
    
    std::ofstream spec((training_path + "/profman.cache").c_str());
    if (!spec.is_open()) {
        return E_NOTFOUND;
    }

    spec << static_cast<unsigned>(correlation_cache.size() + confusion_cache.size()) << std::endl;
    
    for (CorrelationCache::const_iterator i = correlation_cache.begin(); i != correlation_cache.end(); ++i) {
        spec << "CORRELATION " << i->first.first << " " << i->first.second << " " << i->second << std::endl;
        if (!spec) {
            return E_IO;
        }
    }

    for (ConfusionCache::const_iterator i = confusion_cache.begin(); i != confusion_cache.end(); ++i) {
        spec << "CONFUSION " << i->first.first << " " << i->first.second << " " << i->second << std::endl;
        if (!spec) {
            return E_IO;
        }
    }
    
    return 0;
}

static int
ReadCacheSpecFile()
{
    int e;
    std::string training_path;
    e = GetTrainingPath(training_path);
    if (FAILURE(e)) {
        return e;
    }
    
    std::ifstream spec((training_path + "/profman.cache").c_str());
    if (!spec.is_open()) {
        return 0;
    }

    unsigned nspec;
    spec >> nspec;
    
    static const size_t BufSz = 1024;
    char buf[BufSz];
        
    while (nspec--) {
        std::string cache_type;
        ProfileId from;
        ProfileId to;
        std::stringstream filename;
        
        spec >> cache_type >> from >> to;
        std::ws(spec);
        
        do {
            memset(buf, 0, BufSz);
            spec.clear(std::ios_base::goodbit);
            spec.getline(buf, BufSz);
            if (spec.bad() || spec.eof()) {
                break;
            }
            filename << buf;
        } while (spec.fail());
        
        if (spec.bad()) {
            return E_IO;
        }
        
        if (cache_type == "CORRELATION") {
            correlation_cache[std::make_pair(from, to)] = filename.str();
        }
        else if (cache_type == "CONFUSION") {
            confusion_cache[std::make_pair(from, to)] = filename.str();
        }
        else {
            return E_INVALID;
        }
    }
    
    return 0;
}


static int
InitializeCache()
{
    if (CacheInitialized()) {
        return 0;
    }
    
    int e = ReadCacheSpecFile();
    if (FAILURE(e)) {
        return e;
    }
    
    cache_init = true;
    return 0;
}


bool
AreProfileCorrelationsCached(const ProfileId &from, const ProfileId &to)
{
    if (!CacheInitialized()) {
        int e = InitializeCache();
        if (FAILURE(e)) {
            return false;
        }
    }
    
    return correlation_cache.find(std::make_pair(from, to)) != correlation_cache.end();
}


int
ImportProfileCorrelations(const ProfileId &from, const ProfileId &to,
                          std::vector<InterProfileCorrelation> &correlations)
{
    if (!CacheInitialized()) {
        int e = InitializeCache();
        if (FAILURE(e)) {
            return e;
        }
    }
    
    CorrelationCache::const_iterator i = correlation_cache.find(std::make_pair(from, to));
    if (i == correlation_cache.end()) {
        return E_NOTFOUND;
    }
    
    std::ifstream ifs(i->second.c_str());
    if (!ifs.is_open()) {
        return E_NOTFOUND;
    }
    
    ifs >> correlations;
    if (!ifs) {
        return E_IO;
    }

    return 0;
}


int
CacheProfileCorrelations(const ProfileManager &profman, const ProfileId &from, const ProfileId &to,
                         const std::vector<InterProfileCorrelation> &correlations)
{
    int e;
    std::string training_path;
    e = GetTrainingPath(training_path);
    if (FAILURE(e)) {
        return e;
    }
    
    std::stringstream filename;
    filename << "/" << from << "-" << to << ".correlation";
    
    std::string filepath = training_path + filename.str();
    std::ofstream corr(filepath.c_str());
    if (!corr.is_open()) {
        return E_IO;
    }

    corr << correlations;
    if (!corr) {
        return E_IO;
    }
    
    correlation_cache[std::make_pair(from, to)] = filepath;
    
    return WriteCacheSpecFile();
}


bool
IsProfileConfusionCached(const ProfileId &from, const ProfileId &to)
{
    if (!CacheInitialized()) {
        int e = InitializeCache();
        if (FAILURE(e)) {
            return false;
        }
    }
    
    return confusion_cache.find(std::make_pair(from, to)) != confusion_cache.end();
}


int
ImportProfileConfusion(const ProfileId &from, const ProfileId &to,
                       std::vector<InterProfileCorrelation> *confusions, unsigned n)
{
    if (!CacheInitialized()) {
        int e = InitializeCache();
        if (FAILURE(e)) {
            return e;
        }
    }
    
    ConfusionCache::const_iterator i = confusion_cache.find(std::make_pair(from, to));
    if (i == confusion_cache.end()) {
        return E_NOTFOUND;
    }
    
    std::ifstream ifs(i->second.c_str());
    if (!ifs.is_open()) {
        return E_NOTFOUND;
    }
    
    for (unsigned i = 0; i < n; i++) {
        ifs >> confusions[i];
        if (!ifs) {
            return E_IO;
        }
    }
    
    return 0;
}


int
CacheProfileConfusion(const ProfileId &from, const ProfileId &to,
                      const std::vector<InterProfileCorrelation> *confusions, unsigned n)
{
    int e;
    std::string training_path;
    e = GetTrainingPath(training_path);
    if (FAILURE(e)) {
        return e;
    }
    
    std::stringstream filename;
    filename << "/" << from << "-" << to << ".confusion";
    
    std::string filepath = training_path + filename.str();
    std::ofstream conf(filepath.c_str());
    if (!conf.is_open()) {
        return E_IO;
    }

    for (unsigned i = 0; i < n; i++) {
        conf << confusions[i];
        if (!conf) {
            return E_IO;
        }
    }
    
    confusion_cache[std::make_pair(from, to)] = filepath;
    
    return WriteCacheSpecFile();
}


}

