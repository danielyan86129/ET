# ifndef ISO_CONTOUR_H
# define ISO_CONTOUR_H

#include "allTypes.h"
#include "commonUtils.h"

class IsoContourExtractor
{
public:
	IsoContourExtractor() {}
	~IsoContourExtractor() {}
	void marchTriangle(
		const vector<TriPoint>& _vts_of_faces, 
		const vector<float>& _vals, 
		const float _iso_value,
		vector<TriPoint>& _iso_vts, 
		//vector<TriEdge>& _iso_edges, 
		vector<TriPoint>& _vts_for_newFaces, 
		vector<float>& _scalar_for_newFaces_vts/*,
		vector<TriFace>& _newFaces*/
		);
public:
	inline static bool above(float _a, float _iso)
	{
		return _a >= _iso;
	}
	inline static bool below(float _a, float _iso)
	{
		return !above(_a, _iso);
	}
	inline static void new_iso_vert(
		const TriPoint& _u, const TriPoint& _v, float _a, float _b, float _iso, TriPoint& _p
		)
	{
		float l = std::abs(_b - _a);
		if (::almost_zero(l, 0.000000001f))
			_p = _u;
		else
			_p = trimesh::mix(_u, _v, std::abs(_a-_iso) / l);
	}
};

# endif