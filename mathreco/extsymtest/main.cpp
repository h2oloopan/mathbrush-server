#include "MathRecognizer.h"

int
main(int argc, char *argv[])
{
	for (scg::UnicodeIterator i = scg::FirstUnicodeSymbol(); i.unicode_value; i = scg::NextUnicodeSymbol(i)) {
		scg::unicode_char a = i.unicode_value;
	}
	
	return 0;
}