#ifndef CONFUSION_H_
#define CONFUSION_H_


#include "profile.h"
#include "feat.h"

namespace scg
{


// parameters are: 1. number of symbols completed, and 2. total number of symbols;
//typedef void (*status_cb_fn)(unsigned, unsigned);

class ProfileManager;

struct InterProfileConfusion
{
    PrototypeId from_id;
    PrototypeId to_id;
    
    InterProfileConfusion() {}
    InterProfileConfusion(const PrototypeId &from, const PrototypeId &to) : from_id(from), to_id(to) {}    
};

int ComputeProfileConfusion(const ProfileManager &profman, const Profile &from, const Profile &to, std::vector<InterProfileCorrelation> *confusion);


}


#endif

