#ifndef MD5_H_
#define MD5_H_


#include <algorithm>
#include <istream>
#include <ostream>
#include <iostream>


namespace scg
{


struct Md5Hash
{
    const static size_t HashSize; // = 16
    
    unsigned char hash[16];
    
	 Md5Hash()
	 {
	 	std::fill(hash, hash+16, 0);
	 }
	 Md5Hash(const Md5Hash &rhs)
	 {
	 	std::copy(rhs.hash, rhs.hash + 16, hash);
	 }
	 Md5Hash &operator=(const Md5Hash &rhs)
	 {
	 	if (this != &rhs) {
			std::copy(rhs.hash, rhs.hash + 16, hash);
		}
		return *this;
	 }

    bool operator==(const Md5Hash &rhs) const
      { return std::equal(hash, hash + HashSize, rhs.hash); }
    bool operator!=(const Md5Hash &rhs) const
      { return !(*this == rhs); }
    bool operator<(const Md5Hash &rhs) const
      { return std::lexicographical_compare(hash, hash + HashSize, rhs.hash, rhs.hash + HashSize); }
};


void ComputeMd5Hash(Md5Hash &hash, unsigned char *data, size_t sz);
//Md5Hash ComputeMd5Hash(unsigned char *data, size_t sz);


std::istream &operator>>(std::istream &is, scg::Md5Hash &hash);
std::ostream &operator<<(std::ostream &os, const scg::Md5Hash &hash);



}


#endif

