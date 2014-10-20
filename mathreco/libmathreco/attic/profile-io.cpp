#include <map>
#include <vector>
#include <iostream>

#include "error.h"
#include "group.h"
#include "info.h"
#include "ink-io.h"
#include "md5.h"
#include "memory.h"
#include "profile.h"
#include "utils.h"

namespace scg
{

const ProfileId Profile::BadId;

const std::string Profile::FileHeader("PROFILE");


static const int NoMoreSamples = 252;


static int
ReadSampleFromInk(std::istream &is, wchar_t &unicode_char, RawStrokeGroup &strokes)
{
    int dummy_unicode;
    is.read(reinterpret_cast<char *>(&dummy_unicode), sizeof(dummy_unicode));
    if (is.fail()) {
        is.clear(std::ios_base::goodbit);
        return NoMoreSamples;
    }
    
    unicode_char = static_cast<wchar_t>(dummy_unicode);
    
    return ReadMicrosoftStrokeGroupInk(is, strokes);
}

static int
ReadSampleFromText(std::istream &is, wchar_t &unicode_char, RawStrokeGroup &strokes)
{
    unsigned short dummy;
    is >> dummy;
    if (is.eof()) {
        return NoMoreSamples;
    }

    unicode_char = static_cast<wchar_t>(dummy);

    try {
        is >> strokes;
    }
    catch (int e) {
        return e;
    }
    
    if (!is) {
        return E_IO;
    }

    return 0;
}


struct StrokeGroupArray
{
    RawStrokeGroup *groups;
    unsigned count;
    unsigned size;
    
    StrokeGroupArray() : groups(DEBUG_NEW RawStrokeGroup[32]), count(0), size(32) {}
};


std::istream &
operator>>(std::istream &is, scg::SymbolInfo &info)
{
    std::string tmp;
    is >> tmp;
    info.set_name(tmp);
    is >> tmp;
    info.set_symbol_attributes(tmp);
    unsigned short unicode;
    is >> unicode;
    info.set_unicode_char(static_cast<wchar_t>(unicode));
    
    return is;
}


std::istream &
operator>>(std::istream &is, scg::Profile &profile)
{
    int error = 0;
    
    std::istream::pos_type start_of_data = is.tellg();
    std::istream::pos_type end_of_data(-1);

    int i;
    for (i = 0; i < 7; i++) {
        if (is.peek() == scg::Profile::FileHeader[i]) {
            is.get();
        }
        else {
            if (i > 0) {
                is.seekg(-i, std::ios_base::cur);
            }
            break;
        }
    }
    
    int (*SampleReaderFunction)(std::istream &, wchar_t &, scg::RawStrokeGroup &) = 0;
    
    if (i == 7) {
        SampleReaderFunction = &scg::ReadSampleFromText;
    }
    else {
        SampleReaderFunction = &scg::ReadSampleFromInk;
    }
    
    std::map<wchar_t, scg::StrokeGroupArray> samples;
    
    for (;;) {
        wchar_t unicode_char;
        
        scg::RawStrokeGroup strokes;
        error = (*SampleReaderFunction)(is, unicode_char, strokes);
        if (error == scg::NoMoreSamples) {
            break;
        }
        if (FAILURE(error)) {
            throw error;
        }
        
        end_of_data = is.tellg();
        
        scg::StrokeGroupArray &group_array = samples[unicode_char];
        ++group_array.count;
        if (group_array.count == group_array.size) {
            scg::RawStrokeGroup *old_array = group_array.groups;
            group_array.size *= 2;
            group_array.groups = DEBUG_NEW scg::RawStrokeGroup[group_array.size];
            std::copy(old_array, old_array + group_array.count, group_array.groups);
            delete[] old_array;
        }
        group_array.groups[group_array.count - 1] = strokes;
    }
    
	 profile.symbols = new Profile::SymbolContainer();

    for (std::map<wchar_t, scg::StrokeGroupArray>::const_iterator i = samples.begin(); i != samples.end(); ++i) {
        const scg::StrokeGroupArray &symbol_samples = i->second;
        
        scg::ProfileSymbol symbol;
        
        symbol.unicode = i->first;

        error = symbol.FillFromSamples(symbol_samples.groups, symbol_samples.groups + symbol_samples.count);
        if (FAILURE(error)) {
            profile.clear();
            goto cleanup;
        }
        
        profile.symbols->push_back(symbol);
    }
    
cleanup:
    if (!is.bad()) {
        std::istream::pos_type data_size = end_of_data - start_of_data;
        char *data_chunk = DEBUG_NEW char[data_size];
        is.seekg(start_of_data);
        is.read(data_chunk, data_size);
        scg::ComputeMd5Hash(profile.id, reinterpret_cast<unsigned char *>(data_chunk), data_size);
		  unsigned b = 1;
        delete[] data_chunk;
		profile.tag_prototypes();
    }
    else {
		throw E_INVALID;
    }
    
    for (std::map<wchar_t, scg::StrokeGroupArray>::const_iterator i = samples.begin(); i != samples.end(); ++i) {
        delete[] i->second.groups;
    }
    
    if (error) {
        throw error;
    }
    
    return is;
}


std::ostream &
operator<<(std::ostream &os, const scg::Profile &profile)
{
    return os;
}


std::istream &
operator>>(std::istream &is, scg::PrototypeId &id)
{
    std::ws(is);
    is.get(); // skip '<'
    
    is >> id.profile;
    is.get(); // skip ';'
    if (!is) {
        throw E_IO;
    }

    is >> id.entry;
    is.get(); // skip ';'
    if (!is) {
        throw E_IO;
    }

    is >> id.prototype;
    is.get(); //skip '>'
    if (!is) {
        throw E_IO;
    }
    
    return is;
}

std::ostream &
operator<<(std::ostream &os, const scg::PrototypeId &id)
{
    return os << "<" << id.profile << ";" << id.entry << ";" << id.prototype << ">";
}


int
OpenProfileFile(const std::string &path, std::ifstream &ifs)
{
    std::ios_base::openmode mode = std::ios_base::in;

	/*
	std::string training_path;
	GetTrainingPath(training_path);
	training_path += '/';
	std::string file_path = training_path.append(path);
	*/
    ifs.open(path.c_str(), mode);
    if (!ifs.is_open()) {
        return E_NOTFOUND;
    }
    
    std::string dummy;
    ifs >> dummy;
    if (dummy != Profile::FileHeader) {
        mode |= std::ios_base::binary;
    }
    
    ifs.close();
    
    ifs.open(path.c_str(), mode);
    if (!ifs.is_open()) {
        return E_NOTFOUND;
    }
    
    return 0;
}


}

