#ifndef PROFILE_H_
#define PROFILE_H_


#include <fstream>
#include <istream>
#include <string>
#include <vector>

#include "algorithm.h"
#include "debug.h"
#include "group.h"
#include "interp.h"
#include "md5.h"
#include "merge.h"
#include "norm.h"
#include "recodefs.h"
#include "ref.h"
#include "stroke.h"


namespace scg
{



struct PrototypeId;
class Profile;


std::istream &operator>>(std::istream &is, Profile &profile);
std::ostream &operator<<(std::ostream &os, const Profile &profile);

std::istream &operator>>(std::istream &is, PrototypeId &id);
std::ostream &operator<<(std::ostream &os, const PrototypeId &id);


typedef Md5Hash ProfileId;


struct PrototypeId
{
    ProfileId profile;
    unsigned entry;
    unsigned prototype;
    
    bool operator==(const PrototypeId &rhs) const
      { return (prototype == rhs.prototype) && (entry == rhs.entry) && (profile == rhs.profile); }

    bool operator!=(const PrototypeId &rhs) const
      { return !(*this == rhs); }
    
    bool operator<(const PrototypeId &rhs) const
    {
        return (profile < rhs.profile)
            || ((profile == rhs.profile) && ((entry < rhs.entry)
                                         || ((entry == rhs.entry) && (prototype < rhs.prototype))));
    }
};


class ProfileSymbol;


struct PrototypeCorrelation
{
    PrototypeId id;
    double correlation;
    
    PrototypeCorrelation() {}
    PrototypeCorrelation(const PrototypeId &Id, double corr) : id(Id), correlation(corr) {}
};


struct Prototype
{
    PrototypeId id;    
    NormalizedStrokeGroup strokes;
    
    const ProfileSymbol *symbol;
    
    std::vector<PrototypeCorrelation> correlations;
    //std::vector<PrototypeCorrelation> confusions[NumMatchers];
       
    Prototype() {}
    explicit Prototype(const PrototypeId &ID) : id(ID) {}
    
    Prototype &operator=(const Prototype &p)
    {
        id = p.id;
        strokes = p.strokes;
        symbol = p.symbol;
        correlations = p.correlations;
        return *this;
    }
};


class ProfileSymbol
{
public:
    //typedef OwnedArrayType(Prototype) PrototypeContainer;
	 typedef std::vector<Prototype> PrototypeContainer;
    typedef PrototypeContainer::iterator iterator;
    
public:
    typedef PrototypeContainer::const_iterator const_iterator;
    
private:
    void tag_prototypes(const ProfileId &profile_id, unsigned symbol)
    {
        unsigned prototype = 0;
        for (iterator i = begin(); i != end(); ++i) {
            i->symbol = this;
            i->id.profile = profile_id;
            i->id.entry = symbol;
            i->id.prototype = prototype++;
        }
    }

    iterator begin()
      { return prototypes.begin(); }
    iterator end()
      { return prototypes.end(); }
    void clear()
      { prototypes.clear(); }
    
public:
    ProfileSymbol() : unicode(0) {}
    
    explicit ProfileSymbol(wchar_t unicode_char)
      : unicode(unicode_char)
      {}
      
    ProfileSymbol &operator=(const ProfileSymbol &symbol)
    {
        unicode = symbol.unicode;
        prototypes = symbol.prototypes;
        return *this;
    }
   
    
    wchar_t unicode_char() const
      { return unicode; }
      
    bool operator==(wchar_t code) const
      { return unicode == code; }

    
    const_iterator begin() const
      { return prototypes.begin(); }
    const_iterator end() const
      { return prototypes.end(); }
    size_t size() const
      { return prototypes.size(); }
    bool empty() const
      { return prototypes.empty(); }
    
private:
    wchar_t unicode;

    PrototypeContainer prototypes;

private:
    struct PrototypeCluster
    {
        Prototype prototype;
        unsigned count;
    };

    template <typename IteratorType>
    int FillFromSamples(IteratorType first, IteratorType last)
    {        
		  return 0;
		  /*
        if (first == last) {
            return 0;
        }

		  std::vector<PrototypeCluster> clusters;

        PrototypeCluster pf;
        NormalizedStrokeGroup normalized_strokes = normalize(*first);
        pf.prototype.strokes = subdivide(normalized_strokes);
        pf.count = 2;
        
        clusters.push_back(pf);
        ++first;
        
        // cluster all the samples that will merge together nicely
        while (first != last) {
            NormalizedStrokeGroup normalized_strokes = normalize(*first);
            NormalizedStrokeGroup sample = subdivide(normalized_strokes);
            bool merged = false;
            
            for (std::vector<PrototypeCluster>::iterator i = clusters.begin(); i != clusters.end(); ++i) {
                NormalizedStrokeGroup average = merge(sample, i->prototype.strokes, 1.0 / i->count);
                if (!average.empty()) {
                    i->prototype.strokes = average;
                    ++i->count;
                    merged = true;
                }
            }
            
            // if this sample could not be added to a cluster, create a new one
            if (!merged) {
                PrototypeCluster pf;
                pf.prototype.strokes = sample;
                pf.count = 2;
                clusters.push_back(pf);
            }
            
            ++first;
        }
        
        // POSSIBLE TODO: for each cluster, find the sample most closely matching it and select it
        // as the representative sample

        // create the actual prototype array that will fill the symbol based on the results of the above
        // (each cluster representative becomes a prototype)
        for (std::vector<PrototypeCluster>::iterator i = clusters.begin(); i != clusters.end(); ++i) {
            i->prototype.symbol = this;
            prototypes.push_back(i->prototype);
        }
        
        return 0;*/
    }
    
    
private:
    friend class Profile;
    friend class ProfileManager;
    friend std::istream &operator>>(std::istream &, scg::Profile &);
};


class Profile
{
public:
    typedef std::vector<ProfileSymbol> SymbolContainer;
	 typedef ref_t<SymbolContainer, shared_ref_t<owned_release_policy> > SymbolContainerRef;
    
public:
    typedef SymbolContainer::iterator iterator;
    typedef SymbolContainer::const_iterator const_iterator;

    static const ProfileId BadId;
    
    static const std::string FileHeader;
    
    
private:
    int tag_prototypes()
    {
        if (id == BadId) {
            return E_NOTREADY;
        }
        
        unsigned symbol = 0;
        for (iterator i = begin(); i != end(); ++i) {
            i->tag_prototypes(id, symbol++);
        }
        
        return 0;
    }

    void clear()
      { symbols->clear(); }

public:
    Profile()
      {}
    Profile(const Profile &rhs) : id(rhs.id), symbols(rhs.symbols) { }

    Profile &operator=(const Profile &rhs)
    {
	 	if (this != &rhs) {
        id = rhs.id;
        symbols = rhs.symbols;
      } 
        return *this;
    }

    
    bool has_entry(wchar_t code)
      { return scg::find(begin(), end(), code) != end(); }
    
    const ProfileId &Id() const
      { return id; }
     
    
    iterator begin()
      { return symbols->begin(); }
    iterator end()
      { return symbols->end(); }
    const_iterator begin() const
      { return symbols->begin(); }
    const_iterator end() const
      { return symbols->end(); }
    size_t size() const
      { return symbols->size(); }
    bool empty() const
      { return symbols->empty(); }

    bool operator==(const Profile &rhs) const
      { return (id != BadId)
		      && (id == rhs.id); }
    bool operator!=(const Profile &rhs) const
      { return !(*this == rhs); }  
    
public:
    ProfileId id;
    
    SymbolContainerRef symbols;
    
private:
    friend class ProfileManager;
    friend std::istream &operator>>(std::istream &, scg::Profile &);
};


int OpenProfileFile(const std::string &path, std::ifstream &test);


}


#endif

