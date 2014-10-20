#ifndef PROFCACHE_H_
#define PROFCACHE_H_


#include "profman.h"


namespace scg
{


bool AreProfileCorrelationsCached(const ProfileId &from, const ProfileId &to);
int ImportProfileCorrelations(const ProfileId &from, const ProfileId &to,
                              std::vector<InterProfileCorrelation> &correlations);
int CacheProfileCorrelations(const ProfileManager &profman, const ProfileId &from, const ProfileId &to,
                             const std::vector<InterProfileCorrelation> &correlations);

bool IsProfileConfusionCached(const ProfileId &from, const ProfileId &to);
int ImportProfileConfusion(const ProfileId &from, const ProfileId &to,
                           std::vector<InterProfileCorrelation> *confusions, unsigned n);
int CacheProfileConfusion(const ProfileId &from, const ProfileId &to,
                          const std::vector<InterProfileCorrelation> *confusions, unsigned n);


}


#endif

