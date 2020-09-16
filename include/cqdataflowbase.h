#pragma once
///////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2018-2019 CacheQ Systems Inc. All rights reserved.
//
// A base fifo class to access AlphaData designs...
///////////////////////////////////////////////////////////////////////////////
#include <map>

namespace cacheq
{
    class CQDataFlowBase : public FPGAInterface
    {
    protected:
        // The dataflow id...
        static const uint64_t DATA_FLOW_ID = 0x00CAC430DA7AF703;

		// The lowest compatible hardware version.
		static const uint32_t COMPAT_MAJOR = 5;
		static const uint32_t COMPAT_MINOR = 25;

        // Data flow addresses.
        static const uint64_t DF_ID = 0x000000000;  // Data flow ID address
        static const uint64_t DF_VER = 0x000000008; // Data flow design
                                                    // version number (MM.mm)
        static const uint64_t DF_CTRL = 0x00FF0000;    // Data flow control
                                                        //  (reset: set bit 0 to
                                                        //   a '1', then '0'
                                                        //   external module reset: reset bit 1 to
                                                        //   a '1', then '0'
                                                        //   pr reconfig: set bit 2 to
                                                        //   a '1', run PR, then '0')
        static const uint64_t DF_ANAME = 0x000000010;   // Data flow ASCII name. 64 bits wide, only
                                                        //   8 characters
        static const uint64_t DF_INFO	 = 0x000000018;    // Fifo info... DF_INFO_1
        static const uint64_t DF_FIFO_ID = 0x000000020;    // More fifo info... DF_INFO_2

        // Some common masks...
        static const uint64_t BYTE0 = 0x00000000000000FF;
        static const uint64_t BYTE1 = 0x000000000000FF00;
        static const uint8_t  BYTE1_SHIFT = 8;
        static const uint64_t BYTE2 = 0x0000000000FF0000;
        static const uint8_t  BYTE2_SHIFT = 16;
        static const uint64_t BYTE3 = 0x00000000FF000000;
        static const uint8_t  BYTE3_SHIFT = 24;

        static const uint64_t BYTE01 = 0x000000000000FFFF;
        static const uint64_t BYTE23 = 0x00000000FFFF0000;
        static const uint8_t  BYTE23_SHIFT = 16;
        static const uint64_t BYTE45 = 0x0000FFFF00000000;
        static const uint8_t  BYTE45_SHIFT = 32;
        static const uint64_t BYTE67 = 0xFFFF000000000000;
        static const uint8_t  BYTE67_SHIFT = 48;

        // Shared memory pool info...
        typedef std::map<int32_t, std::pair<uint64_t, uint64_t> > PooltoAddrMap;
        typedef PooltoAddrMap::iterator                       PooltoAddrMapIter;

    public:
        static const int32_t FIFO_SIZE = 16;
        CQDataFlowBase(uint16_t id, const std::string &fifoname) :
            _numslotsavail(0), _id(id),
            _name(fifoname), _accessaddr(0),
			_pollingthread (nullptr) { }
		virtual ~CQDataFlowBase() { }

        // Do some static initialization. This will only be called once
        //  for all FPGA fifos of a given type.
        virtual void static_initialize();
        // Do some dynamic initialization. This will be called once
        //  per FPGA fifo.
        virtual void initialize();

        virtual void close_board();

        virtual int32_t getnumslots();
        virtual void reset();
        virtual int32_t size() { return FIFO_SIZE; }
        virtual bool is_infifo() const { return _isinfifo; }
        virtual bool full();
        virtual bool empty();

        virtual int32_t getpoolidfromaddr(uint64_t addr);
        virtual uint64_t getfpgapooladdr(int hostpoolid);
        virtual uint64_t getfpgapoolsize(int hostpoolid);

        static void dumpextrahwfifoids();
#ifdef _DEBUG
        virtual uint64_t getrawfifoinavail();
#endif

    protected:
        bool getfifodata();
        void dumpdfinfo();

        typedef std::map<uint16_t, std::pair<bool, uint8_t>> FifoID2idx;
        typedef FifoID2idx::iterator        FifoID2idxIter;
        static FifoID2idx _fifoid2idx;
        void readfifoids();

        uint16_t _id;
        bool _isinfifo;
        std::string _name;
        uint64_t _accessaddr;

        // Polling thread...
        static void polling_thread();
		static std::mutex *_updatefifolock;
		//std::unique_ptr<std::thread> _pollingthread;
		std::thread *_pollingthread;

    public:
        // This is called either by the polling thread or
        // if single threaded, then by generated main()
        static void updateslotcounts();

    protected:
        static std::vector<CQDataFlowBase *> _fpgafifos;
        static std::mutex _slotcountlock;
        static bool _endthread;
        int32_t _numslotsavail;

		// Dataflow IDs
		static uint32_t _majorver;
		static uint32_t _minorver;
		static std::string _dfname;

        // Fifo data...
        static uint32_t _numinfifos;
        static uint32_t _numoutfifos;
        static uint64_t _fifoinavailoffset;
        static uint64_t _fifooutavailoffset;
        static uint64_t _fifoinidoffset;
        static uint64_t _fifooutidoffset;
        static uint64_t _fifoin0offset;
        static uint64_t _fifoout0offset;
        static uint64_t _dfid;
        static uint32_t _hash;

        static bool _initialized;

    private:
        // Don't allow copy!
        CQDataFlowBase(const CQDataFlowBase&);
        CQDataFlowBase& operator=(const CQDataFlowBase&);
    };
}
