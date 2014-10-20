#include "info.h"

#include <algorithm>


namespace scg
{


std::ostream &
operator<<(std::ostream &os, const scg::SymbolInfo &info)
{
    os << info.name() << std::endl;
    os << info.symbol_attributes() << std::endl;
    os << (unsigned short)info.unicode_char() << std::endl;
    return os;
}

bool
operator==(const scg::SymbolInfo &left, const scg::SymbolInfo &right)
{
    return left.name() == right.name();
}

bool
operator!=(const scg::SymbolInfo &left, const scg::SymbolInfo &right)
{
    return left.name() != right.name();
}


bool
operator==(const scg::SymbolInfo &left, const std::string &right)
{
    return left.name() == right;
}




SymbolInfo::SymbolInfo()
    : symbol_type(SYMBOLTYPE_NORMAL), attr_valid(false), unicode(0)
    {}
	
SymbolInfo::SymbolInfo(const std::string name_)
    : symbol_name(name_), symbol_type(SYMBOLTYPE_NORMAL), attr_valid(false), unicode(0)
	{}
	
	
const std::string &
SymbolInfo::symbol_attributes() const
{
	if (attr_valid) {
	    return attribute_string;
	}
	
    if (symbol_type & SYMBOLTYPE_MATH) {
        attribute_string = "math-symbol";
    }
    
    if (symbol_type & SYMBOLTYPE_STACKED) {
        if (!attribute_string.empty()) {
            attribute_string += ',';
        }
        attribute_string += "stacked";
    }
    if (symbol_type & SYMBOLTYPE_CONTAINER) {
        if (!attribute_string.empty()) {
            attribute_string += ',';
        }
        attribute_string += "container";
    }
    if (symbol_type & SYMBOLTYPE_LFENCE) {
        if (!attribute_string.empty()) {
            attribute_string += ',';
        }
        attribute_string += "lfence";
    }
    if (symbol_type & SYMBOLTYPE_RFENCE) {
        if (!attribute_string.empty()) {
            attribute_string += ',';
        }
        attribute_string += "rfence";
    }

    if (attribute_string.empty()) {
        attribute_string = "none";
    }
    
	attr_valid = true;
	return attribute_string;
	
}
	

void
SymbolInfo::set_symbol_attributes(const std::string &s)
{
	parse_attribute_string(s);
	attr_valid = false;
}
	
	
int
SymbolInfo::get_symbol_type() const
{
    return symbol_type;
}

void
SymbolInfo::setattr_normal_symbol()
{
    symbol_type |= SYMBOLTYPE_NORMAL;
    symbol_type &= ~SYMBOLTYPE_MATH;
    attr_valid = false;
}
    
void
SymbolInfo::setattr_math_symbol()
{
    symbol_type |= SYMBOLTYPE_MATH;
    symbol_type &= ~SYMBOLTYPE_NORMAL;
    attr_valid = false;
}

void
SymbolInfo::setattr_stacked_symbol()
{
    symbol_type |= SYMBOLTYPE_STACKED;
    attr_valid = false;
}


void
SymbolInfo::setattr_container()
{
    symbol_type |= SYMBOLTYPE_CONTAINER;
    attr_valid = false;
}


void
SymbolInfo::setattr_lfence()
{
    symbol_type |= SYMBOLTYPE_LFENCE;
    attr_valid = false;
}


void
SymbolInfo::setattr_rfence()
{
    symbol_type |= SYMBOLTYPE_RFENCE;
    attr_valid = false;
}

void
SymbolInfo::setattr_rangeop()
{
    symbol_type |= SYMBOLTYPE_RANGEOP;
    attr_valid = false;
}

void
SymbolInfo::setattr_relop()
{
    symbol_type |= SYMBOLTYPE_RELOP;
    attr_valid = false;
}


void
SymbolInfo::parse_attribute_string(const std::string &s)
{
    std::string::size_type i = 0;
    
    setattr_normal_symbol();
    
    for (;;) {
        std::string::size_type n = s.find_first_of(',', i);
        
        std::string attr = s.substr(i, n - i);
        if (attr == "math-symbol") {
            setattr_math_symbol();
        }
        else if (attr == "stacked") {
            setattr_stacked_symbol();
        }
        else if (attr == "container") {
            setattr_container();
        }
        else if (attr == "lfence") {
            setattr_lfence();
        }
        else if (attr == "rfence") {
            setattr_rfence();
        }
        else if (attr == "rangeop") {
            setattr_rangeop();
        }
        else if (attr == "relop") {
            setattr_relop();
        }
        if (n == std::string::npos) {
            break;
        }
        
        i = n + 1;
    }
    
    attr_valid = false;
}


}

