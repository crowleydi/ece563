#ifndef __MINIMAT_H__
#define __MINIMAT_H__

class minimat {
	private:
		cplx *m_ptr = 0;
		size_t m_cols = 0;
	
		// hide these so they're not public
		minimat(const minimat&);
		minimat& operator=(const minimat&);

	public:
		minimat(size_t rows, size_t cols)
		{
			m_ptr = new cplx[rows*cols];
			m_cols = cols;
		}

		~minimat()
		{
			delete [] m_ptr;
		}

		cplx& operator()(size_t r, size_t c)
		{
			return m_ptr[r*m_cols + c];
		}

		cplx *memptr()
		{
			return m_ptr;
		}
};


class minivect
{
	private:
		cplx *m_ptr = 0;
		size_t m_size = 0;

		minivect(const minivect&);
		minivect& operator=(const minivect&);

	public:
		minivect(size_t s)
		{
			m_ptr = new cplx[s];
			m_size = s;
		}

		~minivect()
		{
			delete [] m_ptr;
		}

		cplx& operator[](size_t i)
		{
			return m_ptr[i];
		}

		cplx *memptr()
		{
			return m_ptr;
		}
};


#endif // __MINIMAT_H__
