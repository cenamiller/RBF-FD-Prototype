#pragma once
///////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2018-2019 CacheQ Systems Inc. All rights reserved.
//
// A fifo class to communicate with a generic Xilinx board...
///////////////////////////////////////////////////////////////////////////////

#include "fpga/cqxdma/cqxdma.h"
#include "cqdataflowbase.h"
namespace cacheq
{
	class CQXdmaInt : public CQDataFlowBase
	{
		static const uint32_t DEVICE_INDEX = 0;     // Assuming only one board
													//  installed.

		static const uint32_t DMA_ENGINE_INDEX = 0;
	public:
		CQXdmaInt(uint16_t id, const std::string &fifoname) :
			CQDataFlowBase(id, fifoname)
		{ }

		virtual void close_board();
		virtual void readwritefifo(void *pBuf, size_t sz, bool write);
		virtual void syncmemwithboard(bool usedma, uint64_t addr, void *pBuf, uint64_t sz, bool write, const std::string &label = "");

		virtual void static_initialize();

	private:
		// Status of last call into API...
		static CQXDMA_RES   _st;
		// Device
		static CQXdma		_xdma;
	};
}