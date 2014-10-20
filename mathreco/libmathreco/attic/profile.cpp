#include "profile.h"
#include "norm.h"
#include "merge.h"
#include "match.h"
#include "recog.h"
#include "interp.h"
#include "stroke-alg.h"
#include "info.h"
#include "functional.h"
#include "feat.h"
#include "error.h"

#include <algorithm>
#include <functional>
#include <cstring>


namespace scg
{


static unsigned
entry_num_symbols(const ProfileEntry &entry)
{
    return entry.nprotos;
}


unsigned
num_symbols(const ProfileEntry &entry)
{
    return entry.nprotos;
}


unsigned
num_symbols(const Profile &profile)
{
    return std::accumulate(profile.begin(), profile.end(), unsigned(0),
                           binary_unary2_composition(std::plus<unsigned>(), std::ptr_fun(&entry_num_symbols)));
}


void
Prototype::clear()
{
    strokes.clear();
    for (iterator data = begin(); data != end(); data++) {
        data->clear();
    }
    ndata = 0;
}

void
ProfileEntry::clear()
{
    for (iterator prototype = begin(); prototype != end(); ++prototype) {
        prototype->clear();
    }
    nprotos = 0;
}


Profile::Profile() : nentries(0) {}
Profile::Profile(const std::string &name) : profile_name(name), nentries(0) {}


void
Profile::clear()
{
    for (iterator entry = begin(); entry != end(); ++entry) {
        entry->clear();
    }
    nentries = 0;
}

    
unsigned
Profile::size() const
{
    return nentries;
}

bool
Profile::empty() const
{
    return size() == 0;
}


ProfileEntry *
Profile::get_entry(wchar_t unicode)
{
    iterator entry = std::find(begin(), end(), unicode);
    return (entry == end()) ? 0 : entry;
}


const ProfileEntry *
Profile::get_entry(wchar_t unicode) const
{
    const_iterator entry = std::find(begin(), end(), unicode);
    return (entry == end()) ? 0 : entry;
}


const ProfileEntry *
Profile::get_entry(unsigned id) const
{
    const ProfileEntry *entry = entries + (id / MAX_ENTRY_PROTOTYPES);
    
    if (entry >= end()) {
        errval = E_INVALID;
        return 0;
    }
    
    return entry;
}


const Prototype *
Profile::get_prototype(unsigned id) const
{
    const ProfileEntry *entry = get_entry(id);
    
    if (!entry) {
        return 0;
    }
    
    const Prototype *prototype = entry->prototypes + (id % MAX_ENTRY_PROTOTYPES);
    
    if (prototype >= entry->end()) {
        errval = E_INVALID;
        return 0;
    }
    
    return prototype;
}


Profile::const_iterator
Profile::begin() const
{
    return entries;
}


Profile::const_iterator
Profile::end() const
{
    return entries + nentries;
}
  
Profile::iterator
Profile::begin()
{
    return entries;
}


Profile::iterator
Profile::end()
{
    return entries + nentries;
}


unsigned
Profile::add_prototype(const SymbolInfo &info, RawStrokeGroup &strokes)
{
    unsigned id = 0;
        
    ProfileEntry *entry = std::find(begin(), end(), info);
    if (entry == end()) {
        if (size() == MAX_PROFILE_ENTRIES) {
            return E_OUTOFMEM;
        }
        
        id = static_cast<unsigned>(entry - begin()) * MAX_ENTRY_PROTOTYPES;
        Prototype proto(id);
        proto.strokes = normalize(strokes);
        //trim_ends(proto.strokes);
        proto.strokes = subdivide(proto.strokes);
        
        proto.dataset[proto.ndata++] = strokes;

        ProfileEntry &newentry = entries[nentries++];
        newentry.info = info;
        newentry.prototypes[newentry.nprotos++] = proto;
    }
    else {
        NormalizedStrokeGroup norm = normalize(strokes);
        //trim_ends(norm);
        norm = subdivide(norm);
        
        ProfileEntry::iterator i;
        for (i = entry->begin(); i != entry->end(); ++i) {
            Prototype &proto = *i;
            
            if (proto.size() < MAX_PROTOTYPE_SAMPLES) {
                double weight = static_cast<double>(proto.size()) / (proto.size() + 1);
                reorder_strokes(norm, proto.strokes, DEFAULT_PRUNING_THRESHOLD);
                NormalizedStrokeGroup combined = merge(proto.strokes, norm, weight);
                if (num_strokes(combined)) {
                    /*double max_conf = 0.0;
                    NormalizedStrokeGroup best_match;

                    for (ProfileEntry::const_iterator j = entry->begin(); j != entry->end(); ++j) {
                        const Prototype &cmp_proto = *j;
                        
                        NormalizedStrokeGroup cmp = subdivide(normalize(cmp_proto.strokes));
                        Match match;
                        compare_groups(combined, cmp, match);
                        
                        if (match.num_strokes && match.confidence > max_conf) {
                            max_conf = match.confidence;
                            best_match = cmp;
                        }
                    }*/
                    
                    
                    proto.strokes = combined;
                    //proto.strokes = best_match;
                    proto.dataset[proto.ndata++] = strokes;
                    
                    id = proto.id;

                    break;
                }
            }
        }
        
        if (i == entry->end()) {
            if (entry->size() == MAX_ENTRY_PROTOTYPES) {
                return E_OUTOFMEM;
            }
            else {
                id = static_cast<unsigned>(entry - begin()) * MAX_ENTRY_PROTOTYPES + static_cast<unsigned>(i - entry->begin());
                Prototype proto(id);
                
                proto.strokes = norm;
                proto.dataset[proto.ndata++] = strokes;
                entry->prototypes[entry->nprotos++] = proto;
            }
        }
    }
    
    return id;
}

int
Profile::add_entry(const SymbolInfo &info)
{
    if (has_entry(info)) {
        return 0;
    }
    
    if (size() == MAX_PROFILE_ENTRIES) {
        return E_OUTOFMEM;
    }
    
    entries[nentries++].info = info;
    return 0;
}


int
Profile::remove_entry(const SymbolInfo &info)
{
    iterator i = std::find(begin(), end(), info);
    
    if (i != end()) {
        iterator j;
        for (j = i; j != end() - 1; j++) {
            *j = *(j + 1);
        }
        
        j->clear();
        
        nentries--;
    }
    
    return 0;
}


const std::string &
Profile::name() const
{
    return profile_name;
}


int
Profile::compute_feature_correlation()
{
    static const scg::NormalizedStrokeFeatures TIGHT_PRUNING_THRESHOLD(25.0, 25.0, 15.0, 15.0, 30.0, 2.0);
    
    correlations.clear();
        
    for (scg::Profile::const_iterator entry = begin(); entry != end(); ++entry) {
        for (scg::ProfileEntry::const_iterator proto = entry->begin(); proto != entry->end(); ++proto) {
            for (scg::Profile::const_iterator cmp_entry = begin(); cmp_entry != end(); ++cmp_entry) {
                for (scg::ProfileEntry::const_iterator cmp_proto = cmp_entry->begin(); cmp_proto != cmp_entry->end(); ++cmp_proto) {
                    double conf;

                    scg::NormalizedStrokeGroup strokes;
                    if (num_strokes(proto->strokes) > num_strokes(cmp_proto->strokes)) {
                        scg::NormalizedStrokeGroup substrokes = scg::normalize(proto->strokes.strokes(), num_strokes(cmp_proto->strokes));
                        strokes = scg::reorder_strokes(cmp_proto->strokes, substrokes, TIGHT_PRUNING_THRESHOLD, &conf);
                        if (conf > 0.0) {
                            correlations.push_back(correlation_t(proto->id, cmp_proto->id, conf));
                        }
                    }
                }
            }
        }
    }
    
    return 0;
}


}

