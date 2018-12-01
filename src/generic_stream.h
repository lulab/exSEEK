#ifndef __GENERIC_STREAM_H__
#define __GENERIC_STREAM_H__
#include <cstdio>
#include <zlib.h>
#include <string.h>

namespace generic_stream {

template <typename Stream>
inline Stream open(const char* path, const char* mode);
template <> inline FILE* open(const char* path, const char* mode) { return std::fopen(path, mode); }
template <> inline gzFile open(const char* path, const char* mode) { return gzopen(path, mode); }

template <typename Stream>
inline Stream dopen(int fd, const char* mode);
template <> inline FILE* dopen(int fd, const char* mode) { return fdopen(fd, mode); }
template <> inline gzFile dopen(int fd, const char* mode) { return gzdopen(fd, mode); }

template <typename Stream>
inline int putc(Stream stream, int c);
template <> inline int putc(FILE* stream, int c) { return std::fputc(c, stream); }
template <> inline int putc(gzFile stream, int c) { return gzputc(stream, c); }

template <typename Stream>
inline int getc(Stream stream);
template <> inline int getc(FILE* stream) { return std::fgetc(stream); }
template <> inline int getc(gzFile stream) { return gzgetc(stream); }

template <typename Stream>
inline int puts(Stream stream, const char* s);
template <> inline int puts(FILE* stream, const char* s) { return std::fputs(s, stream); }
template <> inline int puts(gzFile stream, const char* s) { return gzputs(stream, s); }

template <typename Stream>
inline char* gets(Stream stream, char* s, int size);
template<> inline char* gets(FILE* stream, char* s, int size) { return std::fgets(s, size, stream); }
template<> inline char* gets(gzFile stream, char* s, int size) { return gzgets(stream, s, size); }

template <typename Stream>
inline size_t write(void* buf, size_t size, size_t nitems, Stream stream);
template<> inline size_t write(void* buf, size_t size, size_t nitems, FILE* stream) { return std::fwrite(buf, size, nitems, stream); }
template<> inline size_t write(void* buf, size_t size, size_t nitems, gzFile stream) { return gzwrite(stream, buf, size*nitems); }

template <typename Stream>
inline size_t read(void* buf, size_t size, size_t nitems, Stream stream);
template<> inline size_t read(void* buf, size_t size, size_t nitems, FILE* stream) { return std::fread(buf, size, nitems, stream); }
template<> inline size_t read(void* buf, size_t size, size_t nitems, gzFile stream) { return gzread(stream, buf, size*nitems); }

template <typename Stream>
inline int setbuffer(Stream stream, size_t size);
template<> inline int setbuffer(FILE* stream, size_t size) { return std::setvbuf(stream, NULL, _IOFBF, size); }
template<> inline int setbuffer(gzFile stream, size_t size) { return gzbuffer(stream, size); }

template <typename Stream>
inline int close(Stream stream);
template<> inline int close(FILE* stream) { return std::fclose(stream); }
template<> inline int close(gzFile stream) { return gzclose(stream); }

template <typename Stream>
inline int eof(Stream stream);
template <> inline int eof(FILE* stream) { return std::feof(stream); }
template <> inline int eof(gzFile stream) { return gzeof(stream); }

template <class Stream>
inline Stream open_file_or_stdin(const char* path, const char* mode)
{
    if(!strcmp(path, "-"))
        return dopen<Stream>(fileno(stdin), mode);
    else
        return open<Stream>(path, mode);
}

template <class Stream>
inline Stream open_file_or_stdout(const char* path, const char* mode)
{
    if(!strcmp(path, "-"))
        return dopen<Stream>(fileno(stdout), mode);
    else
        return open<Stream>(path, mode);
}

} // end of namespace stream

#endif