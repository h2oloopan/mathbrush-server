#ifndef STREAM_DEFS_H_
#define STREAM_DEFS_H_


#define CHECK_ISTREAM_BASIC(stream) \
		if (stream.eof()) { \
			THROW_ERROR(E_EOF, "CHECK_ISTREAM failed with EOF"); \
		} \
		if (!stream) { \
			THROW_ERROR(E_IO, "CHECK_ISTREAM failed with bad stream"); \
		}

#define CHECK_ISTREAM(stream, error) \
		if (stream.eof()) { \
			THROW_ERROR(E_EOF, error); \
		} \
		if (!stream) { \
			THROW_ERROR(E_IO, error); \
		}


#endif
