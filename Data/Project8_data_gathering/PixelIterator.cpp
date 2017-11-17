#include "PixelIterator.h"



PixelIterator::PixelIterator(){
	Init();
}


PixelIterator::~PixelIterator(){
}

bool PixelIterator::GetPixel(int &x, int &y, int &nx, int &ny) {
	int i = ix++;
	//if we ran out of pixel blocks to render
	if (i >= m_numBlocksX * m_numBlocksY) return false;

	/*if (i == m_numBlocksX * m_numBlocksX) {
		int one = 2;
	}*/
	
	//find the indices of the pixel block
	int bx = i % m_numBlocksX;
	int by = i / m_numBlocksX;
	//find the pixel locations for the top-left pixel of the block
	x = bx * m_blockSize;
	y = by * m_blockSize;

	nx = m_blockSize;
	//if we're on the last column
	if (bx == m_numBlocksX - 1) {
		//check if we need to make it not so wide
		int remainder = ( m_numBlocksX * m_blockSize ) % m_width;
		remainder = m_width % m_blockSize;
		if (remainder > 0) {
			nx = remainder;
		}
	}

	ny = m_blockSize;
	//if we're on the last row
	if (by == m_numBlocksY - 1) {
		//check if we need to make it not so tall
		int remainder = ( m_numBlocksY * m_blockSize ) % m_height;
		remainder = m_height % m_blockSize;
		if (remainder > 0) {
			ny = remainder;
		}
	}

	//x = i % m_width;
	//y = i / m_width;
	return true;
}