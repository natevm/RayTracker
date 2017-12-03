#ifndef _PIXEL_ITER_H_INCLUDED_
#define _PIXEL_ITER_H_INCLUDED_

#include <atomic>	//std::atomic

class PixelIterator
{
public:
	PixelIterator();
	~PixelIterator();
	void Init() { ix = 0; m_width = 0; m_height = 0; m_blockSize = 0; m_numBlocksX = 0; m_numBlocksY = 0; }
	
	void setSize(size_t width, size_t height, size_t blockSize) {
		m_width = width;
		m_height = height;
		m_blockSize = blockSize;	//TODO: make sure this isn't bigger than the smallest dimension of the image
		m_numBlocksX = m_width / m_blockSize + ((m_width % m_blockSize == 0)? 0 : 1);
		m_numBlocksY = m_height / m_blockSize + ((m_height % m_blockSize == 0) ? 0 : 1);
	}

	//pass in references for the pixel block's top-left corner, x and y. Also the dimensions of the block
	bool GetPixel(int &x, int &y, int &nx, int &ny);

private:
	std::atomic<int> ix;				//the next pixel block index
	size_t m_width, m_height;			//the width and height of the final render
	size_t m_blockSize;					//the side length of a square block of pixels for a thread to render
	size_t m_numBlocksX, m_numBlocksY;	//the number of pixel blocks along the image's width and height
};


#endif