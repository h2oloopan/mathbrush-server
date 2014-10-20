#include "symbols.h"
#include "grammar-values.h"
#include "normalize.h"
#include "ink-io.h"
#include "stkutils.h"
#include "utils.h"
#include "error.h"
#include "MathRecognizer.h"
#include "mathrecognizer-private.h"

#ifdef USING_TABLETPC
#include "tpc-group.h"
#endif

#include <fstream>
#include <map>
#include <stdio.h>
#ifndef USING_TABLETPC
#include <dirent.h>
#endif
#include <cstring>
#include <cassert>

namespace scg {

enum {
	SYM_NORMAL = 0,
	SYM_STACKED = 1,
	SYM_CONTAINER = 2,
	SYM_SMALL = 4,
	SYM_DOTTED = 8
};

bool symbol::is_normal() const { return attr == SYM_NORMAL; }
bool symbol::is_stacked() const { return attr & SYM_STACKED; }
bool symbol::is_container() const { return attr & SYM_CONTAINER; }
bool symbol::is_small() const { return attr & SYM_SMALL; }
bool symbol::is_dotted() const { return attr & SYM_DOTTED; }

prototype * symbol::firstproto() { return prots; }
prototype * symbol::nextproto(prototype *cur) { return cur->next; }

prototype::prototype() : owner(0), next(0) { }

symbol::symbol() : unicode(0), attr(0), sid(InvalidSemanticId), rclass((1 << AGGREGATE_CLASS) | (1 << SYMBOL_CLASS)), nt(0), prots(0) { }
symbol::~symbol() { clear(); }

void
symbol::clear() {
	while (prots) {
		prototype *nx = prots->next;
		delete prots;
		prots = nx;
	}
}


static void
get_data_path(std::string &path, const char *name, const char *suffix) {
	path += '/';
	path += name;
	path += suffix;
}

static int
get_db_path(const char *name, const char *suffix, std::string &path) {
	int e = GetProfilePath(path);
	if (FAILURE(e)) return e;
	get_data_path(path, name, suffix);
	return 0;
}

static int
get_userdb_path(const char *name, const char *suffix, std::string &path) {
	int e = GetUserProfilePath(path);
	if (FAILURE(e)) return e;
	get_data_path(path, name, suffix);
	return 0;
}

template <typename T>
static int
open_data_file(std::string &path, T &stream) {
	stream.open(path.c_str());
	return stream.is_open() ? E_OK : E_NOTFOUND;
}

static int
open_db_file(const char *name, const char *suffix, std::ifstream &stream) {
	std::string path;
	int e = get_db_path(name, suffix, path);
	if (FAILURE(e)) return e;
	return open_data_file(path, stream);
}

template <typename T>
static int
open_userdb_file(const char *name, const char *suffix, T &stream) {
	std::string path;
	int e = get_userdb_path(name, suffix, path);
	if (FAILURE(e)) return e;
	return open_data_file(path, stream);
}


struct symdb {
	std::vector<symbol> db;
	std::map<unicode_char, symbol *> db_unicode;
	std::map<std::string, symbol *> db_name;
	std::map<SemanticId, symbol *> db_sid;

	~symdb();

	bool empty() const { return db.empty(); }
	void clear();
	int remove(const std::string &name, bool idx);

	void reindex();

	symbol *findsym_name(const std::string &name);
	symbol *findsym_unicode(unicode_char unicode);
	symbol *findsym_sid(SemanticId sid);

	symbol *firstsymbol();
	symbol *nextsymbol(symbol *cur);
};
static symdb maindb;

symdb::~symdb() { clear(); }

void
symdb::clear() {
	db_unicode.clear();
	db_name.clear();
	db_sid.clear();
	for (symbol *s = firstsymbol(); s; s = nextsymbol(s)) {
		s->clear();
	}
	db.clear();
}

void
symdb::reindex() {
	db_unicode.clear();
	db_name.clear();
	db_sid.clear();
	for (std::vector<symbol>::iterator i = db.begin(); i != db.end(); ++i) {
		db_unicode[i->unicode] = &*i;
		db_name[i->name] = &*i;
		db_sid[i->sid] = &*i;
	}
}

int
symdb::remove(const std::string &name, bool idx) {
	symbol *s = findsym_name(name);
	if (!s || s < &db[0] || s >= &db[db.size()]) {
		return E_NOTFOUND;
	}
	db.erase(db.begin() + (s - &db[0]));
	reindex();
	return 0;
}

template <typename T>
static symbol *
findsym(std::map<T, symbol *> &M, const T &key) {
	typename std::map<T, symbol *>::iterator i = M.find(key);
	return (i == M.end()) ? 0 : i->second;
}

symbol *symdb::findsym_name(const std::string &name) { return findsym(db_name, name); }
symbol *symdb::findsym_unicode(unicode_char unicode) { return findsym(db_unicode, unicode); }
symbol *symdb::findsym_sid(SemanticId sid) { return findsym(db_sid, sid); }

symbol *symdb::firstsymbol() { return &db[0]; }
symbol *symdb::nextsymbol(symbol *cur) { return (&db[0] <= cur && cur < &db[db.size()-1]) ? cur+1 : 0; }


static int
import_symbol(std::istream &in, symbol &s) {
	in >> s.name >> s.mathml >> s.latex >> s.sage;
	std::stringstream ss;
	std::string attr;
	for (;;) {
		in >> attr;
		if (in.eof()) {
			return E_EOF;
		}
		if (!in) {
			return E_INVALID;
		}
		if (attr == "normal") {
			s.attr = SYM_NORMAL;
		}
		else if (attr == "stacked") {
			s.attr |= SYM_STACKED;
		}
		else if (attr == "container") {
			s.attr |= SYM_CONTAINER;
		}
		else if (attr == "small") {
			s.attr |= SYM_SMALL;
		}
		else if (attr == "dotted") {
			s.attr |= SYM_DOTTED;
		}
		else {
			rclass_t rc = tag_to_rclass(attr);
			if (rc == ~0) {
				break;
			}
			else {
				s.rclass |= 1 << rc;
			}
		}
	}

	ss << attr;
	ss >> std::hex >> s.unicode;
	if (!ss) {
		return E_INVALID;
	}
	s.sid = mktermsid(s.unicode);
	return 0;
}

static void
prepare_prototype(prototype &prot) {
	prot.nstrokes = normalize(prot.strokes);
}

static int
read_prototype(std::istream &in, prototype &prot) {
	try {
		in >> prot.strokes;
	}
	catch (int e) {
		return e;
	}
	prepare_prototype(prot);
	return 0;
}

static int
verify_file_header(std::istream &in, const std::string &header, unsigned version) {
	size_t i = 0;
	while (i < header.size() && in.peek() == header[i]) {
		in.get();
		++i;
	}
	if (i < header.size()) {
		return E_INVALID;
	}
	
	unsigned fileversion;
	in >> fileversion;
	if (!in || in.eof()) {
		return E_INVALID;
	}
	else if (version != fileversion) {
		return E_OUTOFDATE;
	}
	return E_OK;
}

const static std::string HEADER_PROTS("MBDB");
const static unsigned VERSION_PROTS = 1;

static int
import_prototypes(std::istream &in, symdb &db) {
	int e = verify_file_header(in, HEADER_PROTS, VERSION_PROTS);
	if (FAILURE(e)) {
		return e;
	}
	symbol *s;
	for (;;) {
		unicode_char unicode;
		in >> std::hex >> unicode >> std::dec;
		if (in.eof()) {
			break;
		}
		if (!in) {
			e = E_INVALID;
			break;
		}
		//VERBOSE(*verb_out << "reading prototype for unicode symbol " << std::hex << unicode << std::dec << std::endl);
		s = db.findsym_unicode(unicode);
		if (!s) {
			/*e = E_NOTFOUND;
			VERBOSE(*verb_out << "profile references unknown unicode symbol " << std::hex << unicode << std::dec << std::endl);
			ERR(e, "profile references unknown unicode symbol " << std::hex << unicode << std::dec);
			break;*/
			prototype dummy;
			e = read_prototype(in, dummy);
			if (FAILURE(e)) {
				break;
			}
		}
		else {
			assert(s->unicode == unicode);

			prototype *prot = new prototype;
			prot->owner = s;
			prot->next = s->prots;
			s->prots = prot;
			e = read_prototype(in, *prot);
			if (FAILURE(e)) {
				break;
			}
		}
	}

	return e;
}

static int
export_prototypes(std::ostream &out, symdb &db) {
	out << HEADER_PROTS << VERSION_PROTS << std::endl;
	if (!out) {
		return E_IO;
	}

	int e = E_OK;
	for (symbol *s = db.firstsymbol(); s; s = db.nextsymbol(s)) {
		for (prototype *p = s->firstproto(); p; p = s->nextproto(p)) {
			out << std::hex << s->unicode << std::dec << std::endl;
			if (!out) {
				return E_IO;
			}
			try {
				out << p->strokes << std::endl;
			}
			catch (int e) {
				return e;
			}
		}
	}
	return 0;
}

const static std::string HEADER_DEFS("MBSD");
const static unsigned VERSION_DEFS = 1;

static int
import_symbol_defs(const char *name, symdb &db) {
	std::ifstream in;
	int e = open_db_file(name, ".symdefs", in);
	if (FAILURE(e)) {
		return e;
	}
	e = verify_file_header(in, HEADER_DEFS, VERSION_DEFS);
	if (FAILURE(e)) {
		return e;
	}
	while (!in.eof()) {
		symbol s;
		e = import_symbol(in, s);
		if (e == E_EOF) {
			break;
		}
		else if (FAILURE(e)) {
			db.clear();
			return e;
		}
		db.db.push_back(s);
	}
	db.reindex();
	return 0;
}

int
symdb_init() {
	std::ifstream in;
	int e = import_symbol_defs("base", maindb);
	if (FAILURE(e)) goto fail;
	e = import_symbol_defs("extra", maindb);
	if (FAILURE(e)) goto fail;
	e = open_db_file("base", ".symbols", in);
	if (FAILURE(e)) goto fail;
	e = import_prototypes(in, maindb);
	if (FAILURE(e)) goto fail;
	return 0;

fail:
	VERBOSE(*verb_out << "  failed with code " << e << std::endl);
	maindb.clear();
	return e;
}

int
symdb_shutdown() {
	maindb.clear();
	return 0;
}

int
symdb_reindex() {
	maindb.reindex();
	return 0;
}

int
symdb_removesym(const std::string &name, bool idx) {
	return maindb.remove(name, idx);
}

symbol *symdb_firstsymbol() { return maindb.firstsymbol(); }
symbol *symdb_nextsymbol(symbol *cur) { return maindb.nextsymbol(cur); }

symbol *symdb_findsymbol_name(const std::string &name) { return maindb.findsym_name(name); }
symbol *symdb_findsymbol_unicode(unicode_char unicode) { return maindb.findsym_unicode(unicode); }
symbol *symdb_findsymbol_sid(SemanticId sid) { return maindb.findsym_sid(sid); }

int
symdb_setuserprofile(const char *name) {
	std::ifstream in;
	int e = open_userdb_file(name, ".symbols", in);
	if (FAILURE(e)) return e;
	return import_prototypes(in, maindb);
}



UnicodeIterator::UnicodeIterator() : unicode_value(0), s(0), db(0) { }

static UnicodeIterator EOD_UnicodeIterator;

DLLDECL UnicodeIterator
FirstUnicodeSymbol() {
	SymbolsDB_Id db = CreateEmptySymbolsDB();
	if (!db) {
		return EOD_UnicodeIterator;
	}
	UnicodeIterator it;
	it.db = db;
	symdb *sdb = reinterpret_cast<symdb *>(db);
	it.s = sdb->firstsymbol();
	if (!it.s) {
		CloseSymbolsDB(db);
		return EOD_UnicodeIterator;
	}
	it.unicode_value = it.s->unicode;
	return it;
}


DLLDECL UnicodeIterator
NextUnicodeSymbol(UnicodeIterator i) {
	if (!i.db || !i.s) {
		return EOD_UnicodeIterator;
	}
	UnicodeIterator nx;
	nx.db = i.db;
	symdb *sdb = reinterpret_cast<symdb *>(i.db);
	nx.s = sdb->nextsymbol(i.s);
	if (!nx.s) {
		CloseSymbolsDB(nx.db);
		return EOD_UnicodeIterator;
	}
	nx.unicode_value = nx.s->unicode;
	return nx;
}



DLLDECL SymbolsDB_Id InvalidSymbolsDB = 0;

DLLDECL SymbolsDB_Id
CreateEmptySymbolsDB() {
	symdb *db = new symdb;
	int e = import_symbol_defs("base", *db);
	if (FAILURE(e)) {
		errval = e;
		delete db;
		return InvalidSymbolsDB;
	}
	return reinterpret_cast<SymbolsDB_Id>(db);
}

DLLDECL SymbolsDB_Id
LoadSymbolsDB(const char *name) {
	std::ifstream in;
	int e = open_userdb_file(name, ".symbols", in);
	if (FAILURE(e)) {
		errval = e;
		return InvalidSymbolsDB;
	}
	SymbolsDB_Id id = CreateEmptySymbolsDB();
	if (id == InvalidSymbolsDB) {
		return id;
	}
	symdb *db = reinterpret_cast<symdb *>(id);
	e = import_prototypes(in, *db);
	if (FAILURE(e)) {
		errval = e;
		delete db;
		return InvalidSymbolsDB;
	}
	return id;
}

DLLDECL int
CloseSymbolsDB(SymbolsDB_Id id) {
	symdb *db = reinterpret_cast<symdb *>(id);
	delete db;
	return 0;
}

DLLDECL int
CloseSymbolsDB(SymbolsDB_Id id, const char *name) {
	symdb *db = reinterpret_cast<symdb *>(id);
	if (!db) {
		return E_INVALID;
	}
	std::ofstream out;
	int e = open_userdb_file(name, ".symbols", out);
	if (FAILURE(e)) {
		return e;
	}
	return export_prototypes(out, *db);
}

DLLDECL int
ClearSymbolPrototypes(SymbolsDB_Id id, unicode_char unicode) {
	symdb *db = reinterpret_cast<symdb *>(id);
	if (!db) {
		return E_INVALID;
	}
	symbol *s = db->findsym_unicode(unicode);
	if (!s) {
		return E_NOTFOUND;
	}
	s->clear();
	return 0;
}

DLLDECL int
CountSymbolPrototypes(SymbolsDB_Id id, unicode_char unicode, size_t *count) {
	symdb *db = reinterpret_cast<symdb *>(id);
	if (!db || !count) {
		return E_INVALID;
	}
	symbol *s = db->findsym_unicode(unicode);
	if (!s) {
		return E_NOTFOUND;
	}
	*count = 0;
	for (prototype *p = s->prots; p; p = p->next) {
		++(*count);
	}
	return 0;
}

DLLDECL int
GetSymbolPrototypes(SymbolsDB_Id id, unicode_char unicode, RawStrokeGroup *inks, size_t *count) {
	symdb *db = reinterpret_cast<symdb *>(id);
	if (!db || !count) {
		return E_INVALID;
	}
	symbol *s = db->findsym_unicode(unicode);
	if (!s) {
		return E_NOTFOUND;
	}

	size_t n = 0;
	for (prototype *prot = s->firstproto(); n < *count && prot; prot = s->nextproto(prot), ++n) {
		inks[n] = prot->strokes.copy();
	}

	*count = n;
	return 0;
}

DLLDECL int
ReleaseSymbolPrototypes(RawStrokeGroup *inks, size_t count) {
	for (size_t n = 0; n < count; ++n) {
		inks[n].clear();
	}
	return 0;
}

#ifdef USING_TABLETPC

DLLDECL int
ReleaseSymbolPrototypes(IInkDisp **inks, size_t count) {
	for (size_t n = 0; n < count; ++n) {
		inks[n]->Release();
	}
	return 0;
}


DLLDECL int
GetSymbolPrototypes(SymbolsDB_Id id, unicode_char unicode, IInkDisp **inks, size_t *count) {
	symdb *db = reinterpret_cast<symdb *>(id);
	if (!db || !count) {
		return E_INVALID;
	}
	symbol *s = db->findsym_unicode(unicode);
	if (!s) {
		return E_NOTFOUND;
	}

	size_t n = 0;
	for (prototype *prot = s->firstproto(); n < *count && prot; prot = s->nextproto(prot), ++n) {
		TPC_StrokeGroup msink;
		inks[n] = make_empty_ink_object();
		if (!inks[n]) {
			return get_error().code;
		}
		convert(msink, prot->strokes.strokes, prot->strokes.strokes + prot->strokes.nstrokes, inks[n]);
		msink.add_to_ink(inks[n]);
	}

	*count = n;
	return 0;
}
#endif


DLLDECL int
AddPrototypeToSymbolsDB(SymbolsDB_Id id, unicode_char unicode, long **x, long **y, unsigned *npoints, size_t nstrokes) {
	symdb *db = reinterpret_cast<symdb *>(id);
	if (!db) {
		return E_INVALID;
	}
	symbol *s = db->findsym_unicode(unicode);
	if (!s) {
		return E_NOTFOUND;
	}

	prototype *prot = new prototype;
	RawStroke *strokes = new RawStroke[nstrokes];
	for (size_t i = 0; i < nstrokes; ++i) {
		RawStroke stk(x[i], y[i], 0, npoints[i]);
#ifdef IPAD_RECOGNIZER
		strokes[i] = triple(stk);
#else
		strokes[i] = stk.copy();
#endif
		stk.reset();
	}
	prot->strokes.set_strokes(strokes, nstrokes);
	prot->owner = s;
	prot->next = s->prots;
	s->prots = prot;
	prepare_prototype(*prot);
	return 0;
}

DLLDECL int
AddPrototypeToSymbolsDB(SymbolsDB_Id id, unicode_char unicode, const RawStrokeGroup &stks) {
	symdb *db = reinterpret_cast<symdb *>(id);
	if (!db) {
		return E_INVALID;
	}
	symbol *s = db->findsym_unicode(unicode);
	if (!s) {
		return E_NOTFOUND;
	}

	prototype *prot = new prototype;
#ifdef IPAD_RECOGNIZER
	prot->strokes = triple(stks);
#else
	prot->strokes = stks.copy();
#endif
	prot->owner = s;
	prot->next = s->prots;
	s->prots = prot;
	prepare_prototype(*prot);
	return 0;
}

#ifdef USING_TABLETPC
DLLDECL int
AddPrototypeToSymbolsDB(SymbolsDB_Id id, unicode_char unicode, IInkDisp *ink) {
	RawStrokeGroup scgink;
	convert(scgink, ink);
	return AddPrototypeToSymbolsDB(id, unicode, scgink);
}
#endif

DLLDECL int
SetUserProfile(const char *name) {
	int e = symdb_setuserprofile(name);
	if (FAILURE(e)) return e;
	return RebuildGrammarData();
}

DLLDECL int
RemoveUserProfile(const char *name) {
	std::string path;
	int e = get_userdb_path(name, ".symbols", path);
	if (FAILURE(e)) return e;
	e = remove(path.c_str());
	if (e != 0) {
		return E_NOTFOUND;
	}
	return 0;
}

DLLDECL void
ReleaseAvailableProfiles(char **bufs, size_t count) {
	for (char **p = bufs; p != bufs + count; ++p) {
		delete[] *p;
	}
}

static std::string profile_suffix(".symbols");
static std::wstring profile_wsuffix(L".symbols");

#ifdef WIN32

DLLDECL int
CountAvailableProfiles(size_t *count) {
	if (!count) return E_INVALID;

	std::wstring path;
	int e = GetUserProfilePathW(path);
	if (FAILURE(e)) return e;
	path += L"/*.symbols";
	
	WIN32_FIND_DATA find_data;
	HANDLE sh = ::FindFirstFile(path.c_str(), &find_data);
	if (sh == INVALID_HANDLE_VALUE) {
		if (::GetLastError() == ERROR_FILE_NOT_FOUND) {
			return 0;
		}
		int e = MAKE_API_ERROR(::GetLastError());
		ERR(e, "error while searching for symbols files");
		return e;
	}

	*count = 1;
	
	for (;;) {
		if (!::FindNextFile(sh, &find_data)) {
			if (::GetLastError() == ERROR_NO_MORE_FILES) {
				break;
			}
			int e = MAKE_API_ERROR(::GetLastError());
			ERR(e, "error while searching for symbols files");
			return e;			
		}
		++(*count);
	}
	
	return 0;
}


DLLDECL int
GetAvailableProfiles(char **bufs, size_t count) {
	std::wstring path;
	int e = GetUserProfilePathW(path);
	if (FAILURE(e)) return e;
	path += L"/*.symbols";

	WIN32_FIND_DATA find_data;
	HANDLE sh = ::FindFirstFile(path.c_str(), &find_data);
	if (sh == INVALID_HANDLE_VALUE) {
		if (::GetLastError() == ERROR_FILE_NOT_FOUND) {
			return 0;
		}
		int e = MAKE_API_ERROR(::GetLastError());
		ERR(e, "error while searching for symbols files");
		return e;
	}

	char **s = bufs;
	size_t L = ::wcslen(find_data.cFileName) - profile_wsuffix.size();
	find_data.cFileName[L] = L'\0';
	size_t mblen = ::wcstombs(0, find_data.cFileName, 0);
	*s = new char[mblen+1];
	::wcstombs(*s, find_data.cFileName, mblen+1);
	++s;
	find_data.cFileName[L] = L'.';
	
	for (;;) {
		if (!::FindNextFile(sh, &find_data)) {
			if (::GetLastError() == ERROR_NO_MORE_FILES) {
				break;
			}
			int e = MAKE_API_ERROR(::GetLastError());
			ERR(e, "error while searching for symbols files");
			return e;			
		}

		size_t L = ::wcslen(find_data.cFileName) - profile_wsuffix.size();
		find_data.cFileName[L] = L'\0';
		size_t mblen = ::wcstombs(0, find_data.cFileName, 0);
		*s = new char[mblen+1];
		::wcstombs(*s, find_data.cFileName, mblen+1);
		++s;
		find_data.cFileName[L] = L'.';

		if (--count == 0) {
			break;
		}
	}
	
	return static_cast<int>(s - bufs);
}

#else // ie. non-Windows platforms

DLLDECL int
CountAvailableProfiles(size_t *count) {
	if (!count) return E_INVALID;

	std::string path;
	int e = GetUserProfilePath(path);
	if (FAILURE(e)) return e;

	::DIR *dir = ::opendir(path.c_str());
	
	*count = 0;
	for (dirent *ent = ::readdir(dir); ent; ent = ::readdir(dir)) {
		unsigned L = std::strlen(ent->d_name);
		if (L > profile_suffix.size()) {
			if (std::strncmp(profile_suffix.c_str(), ent->d_name + L - profile_suffix.size(), profile_suffix.size()) == 0) {
				++(*count);
			}
		}
	}

	return 0;
}


DLLDECL int
GetAvailableProfiles(char **bufs, size_t count) {
	std::string path;
	int e = GetUserProfilePath(path);
	if (FAILURE(e)) return e;

	::DIR *dir = ::opendir(path.c_str());
	
	char **s = bufs;
	for (dirent *ent = ::readdir(dir); ent; ent = ::readdir(dir)) {
		unsigned L = std::strlen(ent->d_name);
		if (L > profile_suffix.size()) {
			if (std::strncmp(profile_suffix.c_str(), ent->d_name + L - profile_suffix.size(), profile_suffix.size()) == 0) {
				size_t len = L-profile_suffix.size();
				*s = new char[len+1];
				memcpy(*s, ent->d_name, len);
				(*s)[len] = 0;
				++s;
				if (--count == 0) {
					break;
				}
			}
		}
	}

	return static_cast<int>(s - bufs);
}


#endif


}
