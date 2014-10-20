#ifndef PROFMAN_H_
#define PROFMAN_H_


#include <map>
#include <vector>

#include "array.h"
#include "info.h"
#include "iterator.h"
#include "profile.h"
#include "symtab.h"


namespace scg
{


struct InterProfileCorrelation;
struct InterProfileConfusion;
/*
struct PrototypeCorrelation
{
    PrototypeId id;
    double correlation;
    
    PrototypeCorrelation() {}
    PrototypeCorrelation(const PrototypeId &Id, double corr) : id(Id), correlation(corr) {}
};*/

// Simple wrapper class to let us return references as value types
// (Note: this is a read-only container)
template <typename ContainerType>
class StlContainerWrapper
{
protected:
    typedef ContainerType container_type;
    
public:
    typedef typename container_type::const_iterator const_iterator;

public:
    explicit StlContainerWrapper(const container_type &cont) : container(cont) {}
    StlContainerWrapper(const StlContainerWrapper &wrap) : container(wrap.container) {}
    
    virtual const_iterator begin() const
      { return container.begin(); }
    virtual const_iterator end() const
      { return container.end(); }
    virtual size_t size() const
      { return container.size(); }
    virtual bool empty() const
      { return container.empty(); }

private:
    const container_type &container;
};


class PrototypeCorrelations : public StlContainerWrapper<std::vector<PrototypeCorrelation> >
{
private:
    typedef StlContainerWrapper<std::vector<PrototypeCorrelation> > ParentType;
    
public:
    explicit PrototypeCorrelations(const ParentType::container_type &corr)
      : ParentType(corr)
      {}
};

class ConfusionPrototypes : public StlContainerWrapper<std::vector<const Prototype *> >
{
private:
    typedef StlContainerWrapper<std::vector<const Prototype *> > ParentType;
    
public:
    explicit ConfusionPrototypes(const ParentType::container_type &conf)
      : ParentType(conf)
      {}
};


struct ProfileManagerTypes
{
    typedef std::vector<Profile> ProfileContainer;
    typedef MultiContainerIterator<ProfileContainer, ProfileContainer::iterator, Profile::iterator> iterator;
    typedef MultiContainerIterator<const ProfileContainer, ProfileContainer::const_iterator, Profile::const_iterator> const_iterator;
};


// A SymbolClass represents all the symbols of a specific class, extracted from all loaded profiles.
class SymbolClass
{
public:
    typedef SelectiveIterator<ProfileManagerTypes::const_iterator, wchar_t> const_iterator;

public:
    SymbolClass() {}
    SymbolClass(const_iterator s, const_iterator e) : first(s), last(e) {}
    
    const_iterator begin() const
      { return first; }
    const_iterator end() const
      { return last; }

private:
    const_iterator first;
    const_iterator last;
};


class ProfileManager
{
public:
    typedef ProfileManagerTypes::ProfileContainer ProfileContainer;
    
public:
    typedef ProfileManagerTypes::const_iterator const_iterator;
    typedef ProfileManagerTypes::iterator iterator;
    
public:
    explicit ProfileManager(const SymbolTable &symtab);
    //~ProfileManager();
    
    int AddProfile(Profile &profile);
    //int RemoveProfile(const std::string &name);

    SymbolClass GetSymbolClassByUnicodeChar(wchar_t unicode) const;    
    const Prototype *GetPrototypeById(const PrototypeId &id) const;
    
    const ProfileSymbol *GetSymbolByName(const std::string &name) const;
    
    const SymbolInfo &GetInfoByPrototypeId(const PrototypeId &id) const;
    const SymbolInfo &GetInfoByPrototype(const Prototype &prototype) const;
    const SymbolInfo &GetInfoBySymbol(const ProfileSymbol &symbol) const;
    const SymbolInfo &GetInfoByUnicodeChar(wchar_t ch) const;

    PrototypeCorrelations GetPrototypeCorrelationsFor(const PrototypeId &id) const;
    PrototypeCorrelations GetPrototypeConfusionsFor(const PrototypeId &id, int matcher) const;

public:
    iterator begin();
    iterator end();
    
    const_iterator begin() const;
    const_iterator end() const;
    
    size_t size() const;
    bool empty() const;


private:
    Prototype *GetPrototype(const PrototypeId &id);
    
    int UpdateProfileCorrelations(const Profile &from, const Profile &to);
    int UpdateProfileConfusion(const Profile &from, const Profile &to);

    int ImportPrototypeCorrelations(const std::vector<InterProfileCorrelation> &correlations);
    int ImportPrototypeConfusion(const std::vector<InterProfileCorrelation> *confusions, unsigned n);
    
public:
    /*struct PrototypeData
    {
        const Prototype *prototype;
        //std::vector<const Prototype *> *confusion_symbols;
        std::vector<PrototypeCorrelation> *correlations;
    };*/
    
    std::map<wchar_t, SymbolInfo> symbol_table;
    //std::map<PrototypeId, PrototypeData> prototypes;
    
    ProfileContainer profiles;
};


int CreateDefaultProfileManager(ProfileManager *& profman);


}


#endif

