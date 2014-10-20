#ifndef PROFILE_IO_H_
#define PROFILE_IO_H_


#include <istream>
#include <ostream>

#include "profile.h"


std::istream &operator>>(std::istream &is, scg::Profile &profile);
std::ostream &operator<<(std::ostream &os, const scg::Profile &profile);


#endif

