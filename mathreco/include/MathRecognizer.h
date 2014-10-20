#ifndef SCG_MATHRECOGNIZER_H_
#define SCG_MATHRECOGNIZER_H_


#include "dlldecl.h"

#include "MathRecoTypes.h"

#include <istream>
#include <ostream>
#include <cstddef>

#include "group.h"


namespace scg
{

DLLDECL void DisableSubtreeFiltering();
DLLDECL void DisableRelationFiltering();

DLLDECL void SetTabletResolution(size_t dpi);

DLLDECL void SetProfilePath(const char *path);
DLLDECL void SetUserProfilePath(const char *path);

DLLDECL void SetGrammar(const char *gmr);
DLLDECL void SetWriterPkgid(unsigned id);

DLLDECL int InitializeRecognizer();
DLLDECL void ShutdownRecognizer();

DLLDECL void EnableAutoTraining();
DLLDECL void DisableAutoTraining();
DLLDECL bool IsAutoTrainingEnabled();

#ifdef USING_TABLETPC
DLLDECL MathRecognizer *CreateMathRecognizer(IInkDisp *ink);
// Returns a recognizer object for the given Microsoft ink.
#endif

DLLDECL MathRecognizer *CreateMathRecognizer(std::istream &in);

DLLDECL MathRecognizer *CreateMathRecognizer(const MathSymbol *symbols, size_t n);
DLLDECL MathRecognizer *CreateMathRecognizer(const RawStrokeGroup &strokes);
// Returns a recognizer object for SCG-format ink

DLLDECL MathRecognizer *CreateMathRecognizer();
// Returns a recognizer object containing no ink

DLLDECL const grammar *GetMathGrammar();

// The following functions control debug output.  They are currently global settings
// but will eventually be brought into the MathRecognizer object so that output targets
// can be selected on a per-recognizer basis.


DLLDECL void NoVerboseOutput();
// Instructs the recognizer to not generate debug output.

DLLDECL void VerboseOutputToStream(std::ostream &os);
// Instructs the recognizer to write debug output to the specified standard stream object
// (Note: do not destroy the stream before destroying all recognition objects)

DLLDECL void VerboseOutputToFile(const char *filename);
// Instructs the recognizer to write debug output to a file with the specified path.

DLLDECL void SetVerbosity(int level);
// Choose the amount of debug output generated.  Higher is more.
// This currently has no effect.


DLLDECL ExpressionTree *CreatePlaceholderExpression();
DLLDECL ExpressionTree *CreateBlankExpression();


// The following declarations deal with profile management

typedef ptrdiff_t SymbolsDB_Id;

struct symbol;
struct UnicodeIterator {
	unsigned short unicode_value;
	UnicodeIterator();
private:
	symbol *s;
	SymbolsDB_Id db;
	friend DLLDECL UnicodeIterator FirstUnicodeSymbol();
	friend DLLDECL UnicodeIterator NextUnicodeSymbol(UnicodeIterator);
};


// Functions to iterate over known unicode values. end-of-symbols
// is indicated by a unicode_value of 0.
DLLDECL UnicodeIterator FirstUnicodeSymbol();
DLLDECL UnicodeIterator NextUnicodeSymbol(UnicodeIterator i);


DLLDECL extern SymbolsDB_Id InvalidSymbolsDB;

DLLDECL int CountAvailableProfiles(size_t *count);
DLLDECL int GetAvailableProfiles(char **bufs, size_t count);
DLLDECL void ReleaseAvailableProfiles(char **bufs, size_t count);

DLLDECL int CountSymbolPrototypes(SymbolsDB_Id id, unicode_char symbol, size_t *count);
DLLDECL int GetSymbolPrototypes(SymbolsDB_Id id, unicode_char unicode, RawStrokeGroup *inks, size_t *count);
DLLDECL int ReleaseSymbolPrototypes(RawStrokeGroup *inks, size_t count);

#ifdef USING_TABLETPC
DLLDECL int GetSymbolPrototypes(SymbolsDB_Id db, unicode_char unicode, IInkDisp **inks, size_t *count);
DLLDECL int ReleaseSymbolPrototypes(IInkDisp **inks, size_t count);
#endif

// the next two functions return InvalidSymbolsDB if the db cannot be created
DLLDECL SymbolsDB_Id CreateEmptySymbolsDB();
DLLDECL SymbolsDB_Id LoadSymbolsDB(const char *name);


DLLDECL int CloseSymbolsDB(SymbolsDB_Id id);
DLLDECL int CloseSymbolsDB(SymbolsDB_Id id, const char *profile_name);


DLLDECL int ClearSymbolPrototypes(SymbolsDB_Id id, unicode_char unicode_value);

DLLDECL int AddPrototypeToSymbolsDB(SymbolsDB_Id id, unicode_char unicode, long **x, long **y, unsigned *npoints, size_t nstrokes);
DLLDECL int AddPrototypeToSymbolsDB(SymbolsDB_Id id, unicode_char unicode, const RawStrokeGroup &strokes);
#ifdef USING_TABLETPC
DLLDECL int AddPrototypeToSymbolsDB(SymbolsDB_Id id, unicode_char unicode, IInkDisp *ink);
#endif

DLLDECL int SetUserProfile(const char *profile_name);
DLLDECL int RemoveUserProfile(const char *profile_name);


}


#endif
