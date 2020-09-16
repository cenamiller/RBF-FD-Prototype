#pragma once
// Some compilers (questa, for instance) don't do pragma once
#ifndef ___cq___mem_h___
#define ___cq___mem_h___

///////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2018-2019 CacheQ Systems Inc. All rights reserved.
//
// mem.h -- memory cache controller simulation model
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <stdarg.h>
#include <mutex>
#include <algorithm>
#include <memory>
#include <vector>
#include <map>

namespace cacheq{

typedef int16_t MemoryPoolId;   // probably never have more than 128, but hey?
enum { POOL_UNKNOWN = 0, 
       POOL_CONST = 1,
       POOL_GLOBAL = 2, 
       POOL_LAST = 63 };

class CacheController;

class MemoryPool
{
    friend class MemController;

public:
    ////////////////////////////////////////////////////////////////////////////
    // Do pointer arithmetic
    class PTR {
    public:
        union { 
            uint64_t        i; 
            char*           p;
        };
    public:
        PTR()   { i = 0; p = 0; }       // make sure the whole thing is zero
        PTR(void* _p) { i = 0; p = (char*)_p; }
        PTR(int64_t _i) { i = _i; }
        PTR(uint64_t _i) { i = _i; }
    };

public:
    // Maximum number of pools supported
    enum { MAXNUMPOOLS = 256 };

    // Granularity of memory pool allocation size.  This is typically the size
    // of the blocks that will be transmitted between CEs.
    enum { POOL_ALLOC_GRANULARITY = 1024 };

public:
    ////////////////////////////////////////////////////////////////////////////
    MemoryPool(int id, size_t maxsize);

    ~MemoryPool();
    ////////////////////////////////////////////////////////////////////////////
    // access

    // Note that I'm assuming here that the object here is allocated at the 
    // beginning of the memory pool so that the base of the pool is the address
    // of this object.
    char*       get_base() const        { return (char*)this; }
    int         get_id() const          { return id; }
    char*       get_end() const         { return get_base() + size; }
    char*       get_maxend() const      { return get_base() + max_size;}
    size_t      get_size()              { return size; }

    size_t      get_alloced_size() const{ return endofmallocedmem; }
    void*       get_endofmallocedmem() const 
                                        { return (char*)this+endofmallocedmem; }
    void        set_endofmallocedmem(void* end)
                                        { endofmallocedmem = 
                                             std::max<int64_t>(endofmallocedmem,
                                                      (char*)end-(char*)this); }
protected:
    // Make everything int64 so that alignment is on 64-bit boundaries for
    // platform compatibility.
    // IMPORTANT NOTE:  NOTHING in here can be a pointer because on other CEs 
    //                  the base address may change!!!!
    int64_t     id;                     // Pool id.
    int64_t     size;                   // current allocation 
    int64_t     endofmallocedmem = 0;   // Offset of highest byte malloced.

public:
    void*       sbrk(size_t newsize);   // Extend the size of this pool

protected:
    int64_t     max_size;               // For now, I preallocate the entire
                                        // pool and this is the max size of
                                        // the pool
protected:
    static MemoryPool*  pools[MAXNUMPOOLS];
    static int          numpools;
public:
    static std::mutex   mut;            // the cq_malloc functions use this
    const char* synchronizing = nullptr;// Set during sync process

public:
    static int get_num_pools()          { return numpools; }

    static MemoryPool* get_memory_pool(int poolid);
    static MemoryPool* get_memory_pool_by_addr(void* p);

    ////////////////////////////////////////////////////////////////////////////
    // Canonical address format has region and pool in the upper 16 and 8 bits, 
    // respectively, with the remaining bits as the offest.  These convert 
    // pointers to between local (logical) address and canonical address.
    typedef uint64_t CanonPtr;
    enum {POOLSHIFT = 40};

    static CanonPtr canonicalize_pointer(void* p, int poolid = POOL_UNKNOWN);

    static char* localize_pointer(CanonPtr p, int poolid = POOL_UNKNOWN,
                                  bool createpool = false);

public:
    // This is in mymalloc.cpp to return the correct type of memory pool
    static MemoryPool* get_new_memory_pool(int id, size_t size = 0);
};

///////////////////////////////////////////////////////////////////////////////
// Base CacheController class
class DRAMCacheController;
class QRAMController;
class Opcode;
class MemController
{
public:
    class CCPrivate;

    typedef int64_t     Time;
    typedef int64_t     Addr;

public:
    // Base functions
    Time access(void* address, bool write, Time t, const Opcode* op)
    {
        MemoryPool::PTR ptr;
        ptr.i = (char*)address - pool->get_base();
        return access(ptr.i, write, t, op);
    }

    void set_log_file(FILE* _log)            { logfile = _log; }
                                        
    // given a poolid, return the cach controller
    static inline MemController* getcontroller(int poolid)
    {
        if (poolid >= controllers->size())
            return nullptr;
        return (*controllers)[poolid].get();
    }

    inline void setcontroller(int poolid) 
    {
        while (controllers->size() <= poolid)
            controllers->push_back(nullptr);
        (*controllers)[poolid].reset(this);
    }

    static inline void alloc_controllers_array()
    {
       // initialize controller array
       if (controllers==nullptr) {
          controllers = new std::vector<std::shared_ptr<MemController>>;
       }
    }

    inline static void setportmap(const std::map<int, 
                                                  const std::vector<int>>* pm) 
    {
        portmap = pm; 
    }

    static inline void printmemstats(bool force=false)
    {
        if (force || getenv("CQ_CACHE_STATS") || getenv("CQ_MEM_STATS"))
        for(auto c: *controllers) {
            if (c) {
                c->print_stats();
            }
        }
    }

    void log(const char* fmt, ...) const   // print to log
    {
        if (!logfile)
           return;
        va_list ap;
        va_start(ap, fmt);
        vfprintf(logfile, fmt, ap);
        va_end(ap);
        fflush(logfile);
    }

public:
    // Overloaded functions

    // Compute the time the access returns given an input time.  Note that
    // I haven't implemented separate sync vs. data time.
    virtual Time access(Addr address, bool write, Time t, const Opcode* op) = 0;

    virtual void print_stats() const = 0;

    virtual ~MemController()      { if (controllers) delete controllers; } 

protected:
    MemoryPool*         pool = nullptr;

    static FILE*        logfile;

    // List of all controllers allocated so far.
    static std::vector<std::shared_ptr<MemController>> *controllers;

    static const std::map<int, const std::vector<int>>* portmap;

protected:
    ///////////////////////////////////////////////
    // helper class for Cache and QRAM
    ///////////////////////////////////////////////
    class AccessMap {                   // Need one of these per stripe.
    public:
        enum { TIMEWINDOW = 16384 };    // how far time can get out of sync
        int8_t          maxaccess = 1;  // max of accesses allowed per cycle
                                        // static so I can call default constr
        int8_t          count[TIMEWINDOW];// circular buffer to count
                                        // simultaneous accesses.
        Time            latestt = -1;   // Latest time in the time window.

        // move time to avoid too many simultaneous accesses.
        Time  shiftaccesstime(Time  t);

        bool            debug = false;
        bool            countmultiplecollisions = true;

    public:
    };

};

///////////////////////////////////////////////////////////////////////////////
// DRAM
class DRAMCacheController: public MemController
{
public:
    typedef int64_t     Time;
    typedef int64_t     Addr;

public:
    // Default cache size; can be overridden by env var
    enum { CACHESIZE = 256 * 1024 };

    // Maximum number of simultaneous accesses that can be handled per cache
    // clock cycle
    enum { ACCESSESPERSTRIPE = 2 };     // Number of simultaneous accesses
    enum { STRIPESHIFT = 3 };           // Each stripe is 64-bits wide.
    enum { NSTRIPES = CACHESIZE / 8192 }; // For BRAM-based cache, 8K/stripe.

    // We're gonna use the same line length as the CPUs do.  Shorter lines will
    // probably do better, but it will be application-dependent
    enum { LINELENGTH = 64 };

    // How long does it take to load a line if a cache miss occurs?  Given
    // memory rates of 2+GHz, this should be pessimistic
    enum { LINELOADDELAY = 50 };        // This is in addition to cache delay

    // Saves might be a bit longer, but probably not enough to matter for
    // our simulations.
    enum { LINESAVEDELAY = 50 };        // This is in addition to cache delay

    // Number of associations in set associative.  1 = directmapped
    // Xilinx System Cache does two or four way.
    enum { NOHIT=0, DIRECTMAPPED=1, TWOWAY=2, FOURWAY=4, EIGHTWAY=8 };

public:
    // To keep the implementation invisible, we construct this way
    static MemController* get_new_DRAMCacheController(MemoryPoolId pool, 
                    int cache_size=CACHESIZE,
                    int load_delay=LINELOADDELAY,
                    int save_delay=LINESAVEDELAY,
                    int num_assoc = DIRECTMAPPED,
                    int line_length = LINELENGTH, 
                    int nstripes=NSTRIPES,
                    int max_access = ACCESSESPERSTRIPE);
public:
    // Compute the time the access returns given an input time.  Note that
    // I haven't implemented separate sync vs. data time.
    virtual Time access(Addr address, bool write, Time t, const Opcode* op) = 0;

    virtual void print_stats() const = 0;

    virtual ~DRAMCacheController()      {}
};

///////////////////////////////////////////////////////////////////////////////
// QRAM
class QRAMController: public MemController
{
public:
    // Default QRAM size.  Overridden by GUI
    enum { QRAMSIZE = 16 * 1024 * 1024 }; // Make it overly large
    enum { NSTRIPES = QRAMSIZE / 32768 }; // For URAM-based cache, 32K/stripe.

    // Maximum number of simultaneous accesses that can be handled per cache
    // clock cycle
    enum { ACCESSESPERSTRIPE = 4 };     // Number of simultaneous accesses
    enum { STRIPESHIFT = 3 };           // Each stripe is 64-bits wide.

public:
    static MemController* get_new_QRAMController(MemoryPoolId pool, 
                                           int size = QRAMSIZE,
                                           int max_access = ACCESSESPERSTRIPE);
                   
public:
    // Compute the time the access returns given an input time.  Note that
    // I haven't implemented separate sync vs. data time.
    virtual Time access(Addr address, bool write, Time t, const Opcode* op) = 0;

    virtual void print_stats() const = 0;

    virtual ~QRAMController()      {}
};


} // end of namespace
#endif
