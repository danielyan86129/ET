#include "hybridSkeleton.h"
#include "commonUtils.h"
#include "origGraph.h"
#include "steinerSubdivision.h"

#include <plyall.h>

#include <tuple>
#include <iostream>
#include <vector>
#include <list>
#include <map>

using namespace std;

// <out vert idx, rel order for sorting, poly dual or not>
typedef std::tuple<unsigned, float, bool> out_item;
void print_direct_g(const map< unsigned, list<out_item> >& _directed_graph);
void print_int_list(const vector<unsigned>& _l);

HybridSkeleton::HybridSkeleton()
{

}

HybridSkeleton::~HybridSkeleton()
{

}

void HybridSkeleton::setup(
	const std::shared_ptr<SteinerGraph>& _stg
	)
{
	m_stg = _stg;
	m_stSubdiv = &m_stg->m_stSubdiv;
	//m_stg->m_stEdgeIdx_dual_map.resize(m_stSubdiv->sizeOfStEdges(), -1);// TODO: get rid of this.
	m_orig_g = m_stg->m_origG.get();
}

void HybridSkeleton::reset()
{
	m_stg.reset();

	m_vts.clear();
	m_vts.shrink_to_fit();
	m_edges.clear();
	m_edges.shrink_to_fit();

	//m_edge_id_map.clear();
	
	/*m_poly_faces.clear();
	m_poly_faces.shrink_to_fit();
	m_poly_fromFaceID.clear();
	m_poly_fromFaceID.shrink_to_fit();*/
    m_tri_faces.clear();
	m_tri_faces.shrink_to_fit();
	/*m_tri_fromPolyID.clear();
	m_tri_fromPolyID.shrink_to_fit();*/

	/* information for pruning */
	// values used in pruning for HS elements
	face_bt2bt3.clear();
	face_bt2bt3.shrink_to_fit();
	face_bt2bt3rel.clear();
	face_bt2bt3rel.shrink_to_fit();
	edges_measure.clear();
	edges_measure.shrink_to_fit();
	edges_rel_measure.clear();
	edges_rel_measure.shrink_to_fit();
	
	// ref count and neighboring information for verts and edges
	// model specific ref count
	// doesn't change from prune to prune
	m_ref_vert_const.clear(); 
	m_ref_vert_const.shrink_to_fit();
	m_ref_edge_const.clear();
	m_ref_edge_const.shrink_to_fit();
	m_nb_faces_for_edge.clear();
	m_nb_faces_for_edge.shrink_to_fit();
	m_nb_edges_for_vert.clear();
	m_nb_edges_for_vert.shrink_to_fit();
	// prune specific ref count
	// valid after each prune
	m_ref_vert_per_prune.clear(); 
	m_ref_vert_per_prune.shrink_to_fit();
	m_ref_edge_per_prune.clear();
	m_ref_edge_per_prune.shrink_to_fit();
	// faces to remove flag
	m_to_remove_face.clear();
	m_to_remove_face.shrink_to_fit();
}

void HybridSkeleton::preprocess()
{
	/* 
	** depending on different dual scheme, vts of hybrid skeleton 
	** may consists of different types of vertices:
	** - if dual scheme avoids original triangle vertices,
	**   h.s. vts = orig vts + steiner vts + dual vts
	** - else if the dual scheme is the new one that already includes orig. vts,
	**   h.s. vts = steiner vts + dual vts
	*/ 
	auto sch = m_stg->getBurnScheme();
	const auto& st_subdiv = m_stg->m_stSubdiv;
	if (sch == SteinerGraph::BurnScheme::ORIGINAL_AND_STEINER)
	{
		const auto& all_g_vts = st_subdiv.getAllVts();
		m_vts.reserve(all_g_vts.size() + m_stg->getDualVts().size());
		std::copy(all_g_vts.begin(), all_g_vts.end(), std::back_inserter(m_vts));
		std::copy(m_stg->getDualVts().begin(), m_stg->getDualVts().end(), std::back_inserter(m_vts));
	}
	else if (sch == SteinerGraph::BurnScheme::STEINER_ONLY)
	{
		m_vts.reserve(st_subdiv.sizeOfStVts() + m_stg->getDualVts().size() );
		for (auto i = 0; i < st_subdiv.sizeOfStVts(); ++i)
			m_vts.push_back(st_subdiv.getStVert(i));
		std::copy(m_stg->getDualVts().begin(), m_stg->getDualVts().end(), std::back_inserter(m_vts));
	}
	cout << "skeleton: vts prepared." << endl;

	/*faces are finer triangles from the dual subdivision of each face on MA*/
	m_tri_faces.clear();
	m_tri_faces = m_stg->m_finerTri_byDual;
	cout << "skeleton: faces prepared." << endl;
	
	/*transform fine-tri-vert-id into current vts list*/
	for ( size_t fi = 0; fi < m_tri_faces.size(); ++fi )
	{
		auto f = m_tri_faces[ fi ];
		//cout << fi<<":" << f[ 0 ] << "," << f[ 1 ] << "," << f[ 2 ];
		for ( size_t i = 0; i < 3; ++i )
		{
			f[ i ] = transform_fine_vert_id( f[ i ] );
		}
		m_tri_faces[ fi ] = f;
		//cout << " -> " << f[ 0 ] << "," << f[ 1 ] << "," << f[ 2 ] << endl;
	}

	/*compute area of each face*/
	float a_min = std::numeric_limits<float>::max();
	float a_max = -1.0f;
	m_face_area.resize( m_tri_faces.size() );
	for ( size_t fi = 0; fi < m_tri_faces.size(); ++fi )
	{
		const auto& f = m_tri_faces[ fi ];
		float area = trimesh::len(
			trimesh::trinorm( m_vts[ f[ 0 ] ], m_vts[ f[ 1 ] ], m_vts[ f[ 2 ] ] )
			);
		m_face_area[ fi ] = area;
		a_min = std::min( area, a_min );
		a_max = std::max( area, a_max );
	}
	cout << "skeleton: each face's area computed. [" << a_min << "," << a_max << "]." << endl;

	/*collect edges*/
	m_edges.clear();
	map<TriEdge, int> edge_id_map;
	// part of edges are from existing isolated edges . add them first.
	auto isolate_edges = m_stg->m_origG->getIsolateEdges();
	for ( auto i = 0; i < isolate_edges.size(); ++i )
	{
		edge_id_map[ transform_isolate_edge( isolate_edges[ i ] ) ] = i;
	}
	// the rest are from faces. each face will be converted to edge-rep.
	TriFace f_e_rep;
	for (auto fi = 0; fi < m_tri_faces.size(); ++fi)
	{
		const auto& f = m_tri_faces[fi];
		// convert face to a triple of edge-ids
		for (unsigned i = 0; i < 3; ++i)
		{
			auto e = util::makeEdge(f[i], f[(i+1)%3]);
			auto find_it = edge_id_map.find(e);
			if (find_it != edge_id_map.end())
			{
				f_e_rep[i] = find_it->second;
			}
			else
			{
				auto pre_size = edge_id_map.size();
				edge_id_map[e] = pre_size;
				f_e_rep[i] = pre_size;
			}
		}
		m_tri_faces[ fi ] = f_e_rep;
	}
	m_edges.resize(edge_id_map.size());
	for (auto it = edge_id_map.begin(); it != edge_id_map.end(); ++it)
	{
		const auto& e = it->first;
		m_edges[it->second] = e;
	}
	edge_id_map.clear();
	cout << "skeleton: edges prepared. " << m_edges.size() << endl;

	/* prepare for pruning */
	// init remove tags
	m_removed[1].assign(m_tri_faces.size(), false);
	m_removed[0].assign(m_edges.size(), false);
	m_removed[2].assign(m_vts.size(), false);
	cout << "skeleton: removed-tags init.ed." << endl;

	// build reference count map
	m_ref_vert_const.assign(m_vts.size(), 0);
	m_ref_edge_const.assign(m_edges.size(), 0);
	m_nb_faces_for_edge.assign(m_edges.size(), vector<unsigned>());
	m_nb_edges_for_vert.assign(m_vts.size(), vector<unsigned>());
	// ref count/nb edges for vts
	for (auto e_it = m_edges.begin(); e_it != m_edges.end(); ++e_it)
	{
		//cout << *e_it << endl;
		m_ref_vert_const[(*e_it)[0]] ++;
		m_ref_vert_const[(*e_it)[1]] ++;
		m_nb_edges_for_vert[(*e_it)[0]].push_back(e_it - m_edges.begin());
		m_nb_edges_for_vert[(*e_it)[1]].push_back(e_it - m_edges.begin());
	}
	cout << "skeleton: ref count/nb-edges for vts prepared." << endl;

	// ref count/nb faces for edges
	for (auto fi = 0; fi < m_tri_faces.size(); ++fi)
	{
		const auto& f = m_tri_faces[fi];
		for (auto i = 0; i < 3; ++i)
		{
			m_ref_edge_const[f[i]]++;
			m_nb_faces_for_edge[f[i]].push_back(fi);
		}
	}
	cout << "skeleton: ref count/nb-faces for edges prepared." << endl;
}

//void HybridSkeleton::prepareFaceExtraction(
//	const vector<TriPoint>& dual_vts, 
//	const vector<TriEdge>& dual_edges)
//{
//	//* vts of hybrid skeleton consists of orig vts and dual vts
//	const auto& orig_g_vts = m_orig_g->vts;
//	m_vts.reserve(orig_g_vts.size() + dual_vts.size());
//	std::copy(orig_g_vts.begin(), orig_g_vts.end(), std::back_inserter(m_vts));
//	std::copy(dual_vts.begin(), dual_vts.end(), std::back_inserter(m_vts));
//
//	m_edges.clear();
//	m_edge_id_map.clear();
//	m_edge_fromOrigEdgeID.clear();
//
//	//* edges of hs consists of two parts: 
//	//* one part of edges formed from subdivision
//	//* Note: end of edge should have index offset by orig vts size if the end is a dual vert
//	vector<unsigned> origAndStVts_onE;
//	vector<unsigned> origAndDualVts_onE;
//	m_triEdgeIdx_dualVts_map.clear();
//	m_triEdgeIdx_dualVts_map.resize(m_orig_g->edges.size());
//	for (unsigned ei = 0; ei < m_orig_g->edges.size(); ++ei)
//	{
//		origAndStVts_onE.clear();
//		origAndDualVts_onE.clear();
//		m_stSubdiv->getVtsOnTriEdge(ei, origAndStVts_onE);
//
//		origAndDualVts_onE.push_back(origAndStVts_onE[0]);
//		for (auto st = origAndStVts_onE.begin(); st != origAndStVts_onE.end() - 1; ++st)
//		{
//			auto& cur_ste = util::makeEdge(*st, *(st+1));
//			int stei = m_stSubdiv->getStEdgeIndex(cur_ste);
//			assert(stei >= 0);
//
//			int e_dual = m_stg->m_stEdgeIdx_dual_map[stei];
//			if (e_dual == -1) // no dual vert on this st edge
//				continue;
//
//			origAndDualVts_onE.push_back(e_dual + orig_g_vts.size());
//			m_triEdgeIdx_dualVts_map[ei].push_back(e_dual + orig_g_vts.size());
//		}
//		origAndDualVts_onE.push_back(origAndStVts_onE.back());
//
//		// connect edges between consecutive dual vts, and dual and orig verts
//		for (auto v_it = origAndDualVts_onE.begin(); v_it != origAndDualVts_onE.end() - 1; ++v_it)
//		{
//			auto& pr = make_pair(
//				v_it == origAndDualVts_onE.begin() || v_it == origAndDualVts_onE.end() - 1 ? 
//				true : false,
//				(v_it+1) == origAndDualVts_onE.begin() || v_it == origAndDualVts_onE.end() - 1 ? 
//				true : false
//				);
//			auto& new_e = util::makeEdge(*v_it, *(v_it+1));
//			m_edges.push_back(new_e);
//			m_edge_fromOrigEdgeID.push_back(ei);
//		}
//	}
//	//m_stg->m_stEdgeIdx_dual_map.clear();
//	cout << "subdivided edges from orig edges added to H.S." << endl;
//
//	// another part is dual edges
//	// only keep those edges with poly dual as end whose degree larger than one
//	// this makes the later stage (poly-face extraction) easier
//	// Note: this leaves some poly-dual vertices unreferenced
//	map<unsigned, unsigned> pDual_degree_map;
//	for (auto it = m_stg->m_eDualFaceIdPr_polyDual_map.begin(); it != m_stg->m_eDualFaceIdPr_polyDual_map.end(); ++it)
//	{
//		pDual_degree_map.find(it->second) == pDual_degree_map.end() ?
//			pDual_degree_map[it->second] = 1 : pDual_degree_map[it->second] ++;
//	}
//	vector<decltype(m_stg->m_eDualFaceIdPr_polyDual_map.begin())> iter_to_remove;
//	for (auto it = m_stg->m_eDualFaceIdPr_polyDual_map.begin(); it != m_stg->m_eDualFaceIdPr_polyDual_map.end(); ++it)
//	{
//		auto find_it = pDual_degree_map.find(it->second);
//		if (find_it->second > 1)
//		{
//			//if (it->first.second == )
//			// remember to offset dual vts index ...
//			TriEdge e = util::makeEdge(
//				it->first.first + orig_g_vts.size(), 
//				it->second + orig_g_vts.size()
//				);
//			m_edges.push_back(e);
//			m_edge_fromOrigEdgeID.push_back(-1); // indicating this is a dual edge
//		}
//		else // remove degree-1 poly dual entry from (edge dual, face id -> poly dual) map
//		{
//			iter_to_remove.push_back(it);
//		}
//	}
//	pDual_degree_map.clear();
//	for (auto iter = iter_to_remove.begin(); iter != iter_to_remove.end(); ++iter)
//	{
//		m_stg->m_eDualFaceIdPr_polyDual_map.erase(*iter);
//	}
//	cout << "select dual edges added to H.S." << endl;
//
//	//* build edge index
//	for (unsigned i = 0; i < m_edges.size(); ++i)
//	{
//		/*auto e_to_find = util::makeEdge(491366, 2286810);
//		if (m_edges[i] == e_to_find)
//			cout << "edge "<< e_to_find << " found."<<endl;*/
//		m_edge_id_map[m_edges[i]] = i;
//	}
//}

//void HybridSkeleton::extractPolyFaces()
//{
//	const auto& orig_g = m_orig_g;
//	const auto& st_subdiv = m_stSubdiv;
//	vector<unsigned> vts2traverse;
//	vector<unsigned> face_edge_list;
//	// <order, rel order, visited>
//	map< unsigned, list<out_item> > directed_graph;
//
//	set<unsigned> f_debug;
//	//f_debug.insert(0);
//
//	m_poly_fromFaceID.clear();
//	m_poly_faces.clear();
//	for (unsigned fi = 0; fi < orig_g->faces.size(); ++fi)
//	{
//		//cout << "fi: "<<fi<<endl; //debug
//		const auto& f = orig_g->faces[fi];
//		vts2traverse.clear();
//
//		//* build the neighborship among dual vts & orig vts:
//
//		// 1. scan through the st vts in m_vts_on_face to obtain
//		//    all dual vts together with orig vts on cur face 
//		//    in the same order as m_vts_on_face, call it L
//		for (unsigned i = 0; i < 3; ++i)
//		{
//			unsigned u = f[i];
//			unsigned v = f[(i+1)%3];
//			const auto& e = util::makeEdge(u, v);
//			int ei = orig_g->getEdgeIdx(e);
//			
//			vts2traverse.push_back(u);
//			const auto& cpy = m_triEdgeIdx_dualVts_map[ei];
//			std::copy(cpy.begin(), cpy.end(), std::back_inserter(vts2traverse));
//			if (u > v)
//				std::reverse(vts2traverse.end() - cpy.size(), vts2traverse.end());
//		}
//
//		// 2. now prepare outgoing graph 
//		//    each above vert has outgoing edges reaching other vts.
//		//    relate a relative order to the outgoing vert
//		//    so that they can be ordered
//		unsigned rel_order = 0;
//		directed_graph.clear();
//		set<unsigned> polyDual_set; // debug
//		for (unsigned i = 0; i < vts2traverse.size(); ++i)
//		{
//			unsigned vi = vts2traverse[i];
//				// outgoing: the next vert with rel order i+1
//				directed_graph[vi].push_back(
//					std::make_tuple(vts2traverse[ (i+1)%vts2traverse.size() ], i+1, false));
//			if (vi >= orig_g->vts.size()) // a dual vert
//			{
//				// additional outgoing: (if exists) the poly dual with rel order i+0.5 
//				auto find_it = m_stg->m_eDualFaceIdPr_polyDual_map.find(make_pair(vi - orig_g->vts.size(), fi));
//				if (find_it != m_stg->m_eDualFaceIdPr_polyDual_map.end())
//				{
//					directed_graph[vi].push_back(
//						std::make_tuple(find_it->second + orig_g->vts.size(), i+0.5f, true));
//					// also add outgoing vi to the poly dual
//					directed_graph[find_it->second + orig_g->vts.size()].push_back(
//						std::make_tuple(vi, vts2traverse.size() - i, false));
//
//					polyDual_set.insert(find_it->second + orig_g->vts.size()); // debug
//				}
//			}
//		}
// 		//if (polyDual_set.size() == 1) // debug
// 		//	cout << "face that only has 1 poly dual: " << fi <<endl;
//
//		// 3. sort outgoing vts of the graph entry
//		// according to the relative order given to each 
//		auto compare_outgoing = [] (const out_item& _o1, const out_item& _o2) {
//			return std::get<1>(_o1) < std::get<1>(_o2);
//		};
//		for (auto it = directed_graph.begin(); it != directed_graph.end(); ++it)
//		{
//			it->second.sort(compare_outgoing);
//		}
//
//		// debug: find the degree of the poly dual if it's the only dual of the cur face
//		unsigned cnt_degree = 0;
//		for (auto it = directed_graph.begin(); it != directed_graph.end(); ++it)
//		{
//			if (polyDual_set.size() == 1 && polyDual_set.count(it->first) )
//			{
//				cnt_degree = it->second.size();
//				break;
//			}
//		}
//
//		if (f_debug.count(fi) > 0)
//		{
//			cout << "extracting poly faces for tri face "<<fi<<" when creating HS..." << endl;
//			cout << "vts2traverse: ";
//			print_int_list(vts2traverse);
//			cout << endl;
//			cout << "directed graph: ";
//			print_direct_g(directed_graph);
//			cout << endl;
//		}
//
//		//* trace poly faces using above out-going neighborship structure
//		//  starting from a vert v whose outgoing list is not exhausted,
//		//  trace out a polygon face by following the outgoing vert,
//		//  delete traced vert as tracing continues and 
//		//  record the poly face info,
//		//  delete v if it's been depleted
//		//  go back to the beginning of 2
//
//		/* 
//		for each v in directed_graph:
//		  cur <- v;
//		  do: 
//		    cur_poly.append(cur);
//		    temp <- cur;
//		  
//		    if cur is an orig vert
//		      cur <- cur's first outgoing vert;
//		      just_visited <- cur handle;
//		    else cur is a dual vert
//		      cur <- cur's first outgoing vert after pre in the outgoing list
//		          or the only remaining outgoing vert 
//		      just_visited <- cur handle;
//		    
//		    pre <- temp;
//		    pre's outgoing list.remove(just_visited);
//		  while (next != v)
//		  end do-while
//		  save cur_poly in output poly list;
//		*/
//		unsigned cnt_polyFaces = 0;
//		for (auto vit = vts2traverse.begin(); vit != vts2traverse.end(); ++vit)
//		{
//			unsigned start_v = *vit;
//			unsigned cur = start_v;
//			bool is_poly_dual = false;
//			auto* out_vts = &( directed_graph.find(cur)->second );
//			unsigned pre = cur;
//			if (out_vts->empty())
//				continue;
//			
//			cnt_polyFaces ++;
//			//cout << "tracing from "<<start_v<<endl; // debug
//
//			m_poly_faces.push_back(vector<unsigned>());
//			auto& cur_poly = m_poly_faces.back();
//
//			auto next_v_it = out_vts->begin(); // just init the handle...
//
//			int iter_cnt = 0;
//			int iter_max = 50;
//			do 
//			{
//				cur_poly.push_back(cur);
//				unsigned cur_cpy = cur;
//
//				// trace to the next vert
//				if (!is_poly_dual) // cur is a vert on the face boundary
//				{
//					next_v_it = out_vts->begin();
//					unsigned next_v = std::get<0>(*next_v_it);
//					if (next_v == pre) 
//					{
//						assert(out_vts->size() >= 2);
//						next_v_it++;
//					}
//				}
//				else // cur is a poly dual vert 
//				{
//					// locate pre in the out_vts list
//					auto out_it = out_vts->begin();
//					for (; out_it != out_vts->end(); ++out_it)
//					{
//						if (std::get<0>(*out_it) == pre)
//							break;
//					}
//					
//					// if cannot find pre, then just use the first out vert as the cur
//					// otherwise fetch the out vert after out_it as the cur in the next iteration
//					out_it == out_vts->end() ? 
//						next_v_it = out_vts->begin() 
//						: (++out_it) == out_vts->end() ? 
//						next_v_it = out_vts->begin() : next_v_it = out_it;
//				}
//
//				// prepare for next iteration
//				cur = std::get<0>(*next_v_it);
//				is_poly_dual = std::get<2>(*next_v_it);
//				pre = cur_cpy;
//				out_vts->erase(next_v_it);
//				out_vts = &( directed_graph.find(cur)->second );
//
//				//cout << "cur, pre: "<<cur<<", "<<pre<<endl;
//				iter_cnt ++;
//			} while (cur != start_v /*&& iter_cnt < iter_max*/);
//			
//			// record from-face for above poly faces
//			m_poly_fromFaceID.push_back(fi);
//
//			// print debug face poly faces or if the poly face is too big
//			if (f_debug.count(fi) > 0 /*|| cur_poly.size() > 30*/)
//			{
//				cout << "face "<<fi<<"'s poly face: ";
//				print_int_list(cur_poly);
//				cout << endl;
//			}
//
//			assert(cur_poly.size() >= 3);
//			// turn the poly face into a list of edge indices (instead of vert indices)
//			face_edge_list.clear();
//			face_edge_list.reserve(cur_poly.size());
//			for (unsigned i = 0; i < cur_poly.size(); ++i)
//			{
//				auto e = util::makeEdge(cur_poly[i], cur_poly[ (i+1)%cur_poly.size() ]);
//				auto find_it = m_edge_id_map.find( e );
//
//				if (find_it == m_edge_id_map.end()) // debug
//				{
//					cout << "cannot find edge "<<e <<" when looking up in edge-id-map"<<endl;
//					cout << "the cur poly: ";
//					print_int_list(cur_poly);
//					cout << endl;
//					cout << "from face: "<<fi<<endl;
//					exit(1);
//				}
//				
//				face_edge_list.push_back(find_it->second);
//			}
//			assert(face_edge_list.size() == cur_poly.size());
//			cur_poly = face_edge_list;
//			
//		} // poly faces for cur face extracted
//
//		// assert: if this face has only one poly dual, then 
//		// # poly faces == # degree of the poly dual
//		if (polyDual_set.size() == 1)
//		{	
//			if( cnt_polyFaces != cnt_degree )
//			{
//				cout << fi<<"has only 1 poly dual "<<*polyDual_set.begin()
//					<<" #poly faces "<<cnt_polyFaces<<" #degree "<<cnt_degree<<endl;
//			}
//		}
//		// sanity check: all out-going edges should be depleted
//		for (auto it = directed_graph.begin(); it != directed_graph.end(); ++it)
//		{
//			assert(it->second.empty());
//		}
//	}
//
//	/* prepare for pruning */
//	// init remove tags
//	m_removed[1].assign(m_poly_faces.size(), false);
//	m_removed[0].assign(m_edges.size(), false);
//	m_removed[2].assign(m_vts.size(), false);
//
//	// build reference count map
//	m_ref_vert_per_prune.assign(m_vts.size(), 0);
//	m_ref_edge_per_prune.assign(m_edges.size(), 0);
//	m_nb_faces_for_edge.resize(m_edges.size(), vector<unsigned>());
//	nb_edges_for_vert.resize(m_vts.size(), vector<unsigned>());
//	// ref count/nb edges for vts
//	for (auto e_it = m_edges.begin(); e_it != m_edges.end(); ++e_it)
//	{
//		m_ref_vert_per_prune[(*e_it)[0]] ++;
//		m_ref_vert_per_prune[(*e_it)[1]] ++;
//		nb_edges_for_vert[(*e_it)[0]].push_back(e_it - m_edges.begin());
//		nb_edges_for_vert[(*e_it)[1]].push_back(e_it - m_edges.begin());
//	}
//	// ref count/nb faces for edges
//	for (auto fit = m_poly_faces.begin(); fit != m_poly_faces.end(); ++fit)
//	{
//		for (auto ei_it = fit->begin(); ei_it != fit->end(); ++ei_it)
//		{
//			m_ref_edge_per_prune[*ei_it] ++;
//			m_nb_faces_for_edge[*ei_it].push_back(fit - m_poly_faces.begin());
//		}
//	}
//}

const vector<TriPoint>& HybridSkeleton::getVts() const
{
	return m_vts;
}

const vector<TriEdge>& HybridSkeleton::getEdges() const
{
	return m_edges;
}

const vector<TriFace>& HybridSkeleton::getTriFaces() const
{
	return m_tri_faces;
}

const vector<bool>& HybridSkeleton::getEdgeRemovedTags() const
{
	return m_removed[0];
}

const vector<bool>& HybridSkeleton::getFaceRemovedTags() const
{
	return m_removed[1];
}

//const vector<unsigned>& HybridSkeleton::getTriFromPolyID() const
//{
//	return m_tri_fromPolyID;
//}

void HybridSkeleton::getRemainedVts(vector<int>& _remain_vts_ids) const
{
	// There are some vertices that are unreferenced (e.g. poly-duals
	// that only have degree-1). We ignore those vertices.
	// Therefore the set of vertices of interest here are those vertices
	// referred to by at least one edge originally
	_remain_vts_ids.clear();
	for (size_t i = 0; i < m_vts.size(); ++i)
	{
		if (m_ref_vert_per_prune[i] == 0)
			continue;
		if (!m_removed[2][i])
		{
			_remain_vts_ids.push_back(i);
			if (m_ref_vert_per_prune[i] <= 0)
			{
				cout << "Attention: vert "<<i<<" is not removed but has ref-cnt "<<m_ref_vert_per_prune[i]<<endl;
			}
		}
		else
		{
			if (m_ref_vert_per_prune[i] > 0)
			{
				cout << "Logic Error: vert "<<i<<" is removed but still has ref-cnt "<<m_ref_vert_per_prune[i]<<endl;
			}
		}
	}
}

void HybridSkeleton::getRemainedEdges(vector<TriEdge>& _edges) const
{
	for (unsigned i = 0; i < m_edges.size(); ++i)
	{
		/*if ( !removed[0][i] )
		{
			_edges.push_back(m_edges[i]);
		}*/

		// how about just returning remaining isolated edges, 
		// i.e. those w/o being used by any face
		if ( !m_removed[0][i] && this->m_ref_edge_per_prune[i] == 0)
		{
			_edges.push_back(m_edges[i]);
		}

		// debug: return ref-cnt == 1 edges
		/*if ( m_ref_edge_per_prune[ i ] == 1 )
		{
			_edges.push_back( m_edges[ i ] );
		}*/
		//_edges.push_back( m_edges[ i ] );
	}
}

void HybridSkeleton::getRemainedFaces(vector<TriFace>& _tris) const
{
	for (unsigned fi = 0; fi < m_tri_faces.size(); ++fi)
	{
		if ( !m_removed[1][fi] )
		{
			int vts_of_tri[ 3 ];
			this->to_vts_list( m_tri_faces[ fi ], vts_of_tri );
			_tris.push_back( TriFace( vts_of_tri ) );
		}
	}
}

void HybridSkeleton::getDualEdges( vector<TriEdge>& _edges ) const
{
	_edges.clear();
	_edges.reserve( m_edges.size() );
	for ( auto i = 0; i < m_edges.size(); ++i )
	{
		const auto& e = m_edges[ i ];
		if ( !part_of_orig_edge( i ) )
			_edges.push_back(e);
	}
}

void HybridSkeleton::getEdgeDiffValColor(float* edge_color)
{
	// obtain orig edge value scale
	float max_val_orig = numeric_limits<float>::min(), min_val_orig = numeric_limits<float>::max();
	float max_val_dual = max_val_orig, min_val_dual = min_val_orig;
	for (unsigned i = 0; i < m_edges.size(); ++i)
	{
		if (part_of_orig_edge(i))
		{
			min_val_orig = std::min(edges_measure[i], min_val_orig);
			max_val_orig = std::max(edges_measure[i], max_val_orig);
		}
		else
		{
			min_val_dual = std::min(edges_measure[i], min_val_dual);
			max_val_dual = std::max(edges_measure[i], max_val_dual);
		}
	}

	unsigned k = 0;
	for (unsigned i = 0; i < m_edges.size(); ++i)
	{
		if (m_removed[0][i])
			continue;

		TriColor c;
		if (part_of_orig_edge(i))
			c = util::GetColour(edges_measure[i], min_val_orig, max_val_orig);
		else
			c = util::GetColour(edges_measure[i], min_val_dual, max_val_dual);

		edge_color[2*3*k + 0] = c[0];
		edge_color[2*3*k + 1] = c[1];
		edge_color[2*3*k + 2] = c[2];
		edge_color[2*3*k + 3] = c[0];
		edge_color[2*3*k + 4] = c[1];
		edge_color[2*3*k + 5] = c[2];

		k++;
	}
}

void HybridSkeleton::getEdgeRelDiffValColor(float* edge_color)
{
	// obtain orig edge value scale
	float max_val_orig = -1.0f, min_val_orig = numeric_limits<float>::max();
	float max_val_dual = -1.0f, min_val_dual = min_val_orig;
	for (unsigned i = 0; i < m_edges.size(); ++i)
	{
		if (part_of_orig_edge(i))
		{
			min_val_orig = std::min(edges_rel_measure[i], min_val_orig);
			max_val_orig = std::max(edges_rel_measure[i], max_val_orig);
		}
		else
		{
			min_val_dual = std::min(edges_rel_measure[i], min_val_dual);
			max_val_dual = std::max(edges_rel_measure[i], max_val_dual);
		}
	}

	unsigned k = 0;
	for (unsigned i = 0; i < m_edges.size(); ++i)
	{
		if (m_removed[0][i])
			continue;

		TriColor c;
		if (part_of_orig_edge(i))
			c = util::GetColour(edges_rel_measure[i], min_val_orig, max_val_orig);
		else
			c = util::GetColour(edges_rel_measure[i], min_val_dual, max_val_dual);

		edge_color[2*3*k + 0] = c[0];
		edge_color[2*3*k + 1] = c[1];
		edge_color[2*3*k + 2] = c[2];
		edge_color[2*3*k + 3] = c[0];
		edge_color[2*3*k + 4] = c[1];
		edge_color[2*3*k + 5] = c[2];

		k++;
	}
}

void HybridSkeleton::getRemainedEdgesColor( vector<TriColor>& _e_colors ) const
{
	// obtain orig edge value scale
	float max_val_orig = numeric_limits<float>::min(), min_val_orig = numeric_limits<float>::max();
	float max_val_dual = max_val_orig, min_val_dual = min_val_orig;
	for ( unsigned i = 0; i < m_edges.size(); ++i )
	{
		/*if ( part_of_orig_edge( i ) )
		{
			min_val_orig = std::min( edges_measure[ i ], min_val_orig );
			max_val_orig = std::max( edges_measure[ i ], max_val_orig );
		}
		else*/
		{
			min_val_dual = std::min( edges_measure[ i ], min_val_dual );
			max_val_dual = std::max( edges_measure[ i ], max_val_dual );
		}
	}

	_e_colors.reserve( m_edges.size() );
	unsigned k = 0;
	for ( unsigned i = 0; i < m_edges.size(); ++i )
	{
		if ( m_removed[ 0 ][ i ] )
			continue;

		TriColor c;
		/*if ( part_of_orig_edge( i ) )
			c = util::GetColour( edges_measure[ i ], min_val_orig, max_val_orig );
		else*/
			c = util::GetColour( edges_measure[ i ], min_val_dual, max_val_dual );

		_e_colors.push_back( c );
	}
	_e_colors.shrink_to_fit();
}

void HybridSkeleton::getTriColor(float* tri_color)
{
	float max_val_orig = -1.0f, min_val_orig = numeric_limits<float>::max();

	for (unsigned i = 0; i < m_tri_faces.size(); ++i)
	{
		//poly
	}
}

void HybridSkeleton::assignElementValues(
	const vector<vector<float>>& _orig_persheet_bt2bt3, 
	const vector<vector<float>>& _orig_persheet_bt2bt3rel, 
	const vector<float>& _orig_vts_bt2bt3,
	const vector<float>& _orig_vts_bt2bt3rel,
	const vector<float>& _dual_vts_bt2bt3,
	const vector<float>& _dual_vts_bt2bt3rel,
	const vector<float>& _dual_vts_bt1bt2, 
	const vector<float>& _dual_vts_bt1bt2rel
	)
{
	min_bt2bt3 = numeric_limits<float>::max(); max_bt2bt3 = numeric_limits<float>::min();
	min_bt2bt3rel = min_bt2bt3; max_bt2bt3rel = max_bt2bt3;
	min_bt1bt2 = min_bt2bt3; max_bt1bt2 = max_bt2bt3;
	min_bt1bt2rel = min_bt2bt3; max_bt1bt2rel = max_bt2bt3;

	// derive values for st edges, and dual edges
	edges_measure.resize(m_edges.size());
	edges_rel_measure.resize(m_edges.size());
	float e_value ;
	unsigned num_orig_vts = m_orig_g->vts.size();
	set<unsigned> debug_e;
	//debug_e.insert( 6 );
	for (unsigned ei = 0; ei < m_edges.size(); ++ei)
	{
		// debug
		//cout << "ei: " << ei << endl;
		
		const auto& e = m_edges[ei];
		//int from_edge = m_edge_fromOrigEdgeID[ei];
		int v0, v1;
		float val_0, val_1;
		bool is_dual_0 = m_stg->isDualVertInFineTri( transform_skel_vert_id( e[ 0 ] ), v0 );
		bool is_dual_1 = m_stg->isDualVertInFineTri( transform_skel_vert_id( e[ 1 ] ), v1 );
		if ( part_of_orig_edge(ei) ) // e is part of orig. tri edge
		{
			val_0 = is_dual_0 ? _dual_vts_bt2bt3[v0] : _orig_vts_bt2bt3[v0];
			val_1 = is_dual_1 ? _dual_vts_bt2bt3[v1] : _orig_vts_bt2bt3[v1];
			e_value = std::min(val_0, val_1);
			/*if ( debug_e.count( ei ) )
				cout << "e-value bt2bt3 done." << e_value << endl;*/
			edges_measure[ei] = e_value;
			min_bt2bt3 = std::min(min_bt2bt3, e_value);
			max_bt2bt3 = std::max(max_bt2bt3, max(val_0, val_1));
			/*if ( debug_e.count( ei ) )
				cout << "range of bt2bt3 done." << endl;*/

			val_0 = is_dual_0 ? _dual_vts_bt2bt3rel[v0] : _orig_vts_bt2bt3rel[v0];
			val_1 = is_dual_1 ? _dual_vts_bt2bt3rel[v1] : _orig_vts_bt2bt3rel[v1];
			e_value = std::min(val_0, val_1);
			/*if ( debug_e.count( ei ) )
				cout << "e-value bt2bt3 rel. done." << e_value << endl;*/
			edges_rel_measure[ei] = e_value;
			min_bt2bt3rel = std::min(min_bt2bt3rel, e_value);
			max_bt2bt3rel = std::max(max_bt2bt3rel, e_value);
			/*if ( debug_e.count( ei ) )
				cout << "range of bt2bt3 rel. done." << endl;*/
		}
		else // e is a dual edge
		{
			assert( is_dual_0 && is_dual_1 );
			
			e_value = std::min(_dual_vts_bt1bt2[v0], _dual_vts_bt1bt2[v1]);
			edges_measure[ei] = e_value;
			if ( ::is_valid(e_value) && e_value < numeric_limits<float>::max() )
			{
				min_bt1bt2 = std::min(min_bt1bt2, e_value);
				max_bt1bt2 = std::max(max_bt1bt2, e_value);
			}

			e_value = std::min(_dual_vts_bt1bt2rel[v0], _dual_vts_bt1bt2rel[v1]);
			edges_rel_measure[ei] = e_value;
			if ( ::is_valid(e_value) && e_value < numeric_limits<float>::max() )
			{
				min_bt1bt2rel = std::min(min_bt1bt2rel, e_value);
				max_bt1bt2rel = std::max(max_bt1bt2rel, e_value);
			}
		}
	}
	std::cout << "HS edges' values derived." << endl;

	// derive values for each tri face
	face_bt2bt3.resize(m_tri_faces.size());
	face_bt2bt3rel.resize(m_tri_faces.size());
	for (unsigned i = 0; i < m_tri_faces.size(); ++i)
	{
		// convert cur face from edge-rep to vert-rep
		vector<unsigned> vts_cur_face;
		const auto& f = m_tri_faces[i];
		auto e = m_edges[f[0]];
		vts_cur_face.push_back(e[0]);
		vts_cur_face.push_back(e[1]);
		e = m_edges[f[1]];
		vts_cur_face.push_back(e[0] != vts_cur_face.back() ? e[0] : e[1]);

		auto from_fi = m_stg->whichOrigFaceIsDualFineFaceFrom(i);
		float scalar_bt2bt3 = 0.0f;
		float scalar_bt2bt3rel = 0.0f;
		int cnt_scalar = 0;
		for (unsigned j = 0; j < vts_cur_face.size(); ++j)
		{
			int v = transform_skel_vert_id( vts_cur_face[ j ] );
			if ( m_stg->isDualVertInFineTri(v, v) ) // v is a dual vert
			{
				if ( m_stg->m_is_face_dual[v] >= 0 )
				{
					scalar_bt2bt3 += _dual_vts_bt2bt3[v];
					scalar_bt2bt3rel += _dual_vts_bt1bt2rel[v];
					cnt_scalar ++;
				}
			}
			else // v is a st vert
			{
				int res = m_stg->mapTopo(v, from_fi);
				assert(res >= 0);

				auto tei = (TopoGraph::EdgeIdx)res;
				scalar_bt2bt3 += _orig_persheet_bt2bt3[v][tei];
				scalar_bt2bt3rel += _orig_persheet_bt2bt3rel[v][tei];
				cnt_scalar ++;
			}
		}
		scalar_bt2bt3 /= cnt_scalar;
		scalar_bt2bt3rel /= cnt_scalar;

		face_bt2bt3[i] = scalar_bt2bt3;
		face_bt2bt3rel[i] = scalar_bt2bt3rel;

		min_bt2bt3 = std::min(face_bt2bt3[i], min_bt2bt3);
		max_bt2bt3 = std::max(face_bt2bt3[i], max_bt2bt3);
		min_bt2bt3rel = std::min(face_bt2bt3rel[i], min_bt2bt3rel);
		max_bt2bt3rel = std::max(face_bt2bt3rel[i], max_bt2bt3rel);
	}
	std::cout << "HS faces' values derived." << endl;

	std::cout << "bt2bt3 diff range: ["<<min_bt2bt3<<","<<max_bt2bt3<<"]"<<endl;
	std::cout << "bt2bt3 rel diff range: ["<<min_bt2bt3rel<<","<<max_bt2bt3rel<<"]"<<endl;
	std::cout << "bt1bt2 diff range: ["<<min_bt1bt2<<","<<max_bt1bt2<<"]"<<endl;
	std::cout << "bt1bt2 rel diff range: ["<<min_bt1bt2rel<<","<<max_bt1bt2rel<<"]"<<endl;
}

void HybridSkeleton::setPolyFaceDegenerateThresh(float _t)
{
	m_face_degen_thresh = _t;
}

void HybridSkeleton::setComponentFaceNumberThresh(int _t)
{
	m_cmpnt_num_faces_thresh = _t;
}

void HybridSkeleton::prune( 
	float _f_diff_r, float _f_reldiff_r, float _l_diff_r, float _l_reldiff_r,
	bool _use_inputs_directly,
	bool _remove_small_components
	)
{
	auto sanity_check = [&]()
	{
		for ( unsigned ei = 0; ei < m_edges.size(); ++ei )
		{
			if ( this->m_ref_edge_per_prune[ ei ] == 0 &&
				!m_removed[ 0 ][ ei ] )
			{
				const auto& e = m_edges[ ei ];
				/*cout << "0-ref edge " << ei << e << " not removed! "
					<< ( part_of_orig_edge( ei ) ? "st edge." : "dual edge." ) <<" "
					<< edges_measure[ei] <<". "
					<< "ref of vert at end: "
					<< m_ref_vert_per_prune[ e[ 0 ] ] << "," << m_ref_vert_per_prune[ e[ 1 ] ]
					<< endl;*/
				//break;
			}
		}
	};

	// get various thresholds 
	float bt2bt3_t, bt2bt3rel_t;
	float bt1bt2_t, bt1bt2rel_t;
	if (_use_inputs_directly)
	{
		bt2bt3_t = _f_diff_r;
		bt2bt3rel_t = _f_reldiff_r;
		bt1bt2_t = _l_diff_r;
		bt1bt2rel_t = _l_reldiff_r;
	}
	else
	{
		ratioToAbs(
			_f_diff_r, 
			&bt2bt3_t, 
			_f_reldiff_r, &bt2bt3rel_t, 
			_l_diff_r, &bt1bt2_t, 
			_l_reldiff_r, &bt1bt2rel_t);
	}

	cout << "use input threshold? " << ( _use_inputs_directly ? "true" : "false" ) << endl;
	cout << "bt2bt3_t, bt2bt3rel_t, bt1bt2_t, bt1bt2rel_t: "
		<< bt2bt3_t << ", " << bt2bt3rel_t << ", " << bt1bt2_t << ", " << bt1bt2rel_t << endl;

	set<unsigned> vts_to_debug;
	//vts_to_debug.insert(3807);
	queue<simple_pair> q;

	/*
	01. build reference count map
	02. push all simple pairs in q that are blow the thresholds
	03. iteratively retracting
	*/
	/*01. reset remove tags & ref count for verts/edges */
	m_removed[0].assign(m_edges.size(), false);
	m_removed[1].assign(m_tri_faces.size(), false);
	m_removed[2].assign(m_vts.size(), false);
	m_ref_vert_per_prune = m_ref_vert_const;
	m_ref_edge_per_prune = m_ref_edge_const;
	// for now no faces will be removed
	m_to_remove_face.assign( m_tri_faces.size(), false);

	/* 02. init q: find all simple pairs whose simple cell has a value below threshold */
	cout << "init.ing q ..." << endl;
	for (unsigned ei = 0; ei < m_edges.size(); ++ei)
	{
		if ( m_ref_edge_per_prune[ ei ] == 1 )
		{
			unsigned nb_f = m_nb_faces_for_edge[ ei ][ 0 ];
			/*float val_bt2bt3 = edges_diff[ei];
			float val_bt2bt3rel = edges_reldiff[ei];*/
			if ( face_edge_below_threshold( nb_f, ei, bt2bt3_t, bt2bt3rel_t, bt1bt2_t, bt1bt2rel_t ) )
			{
				q.push( simple_pair( 1, nb_f, ei ) );
			}
		}
	}
	for ( unsigned vi = 0; vi < m_vts.size(); ++vi )
	{
		if ( m_ref_vert_per_prune[ vi ] == 1 )
		{
			unsigned nb_e = m_nb_edges_for_vert[ vi ][ 0 ];
			if ( edge_below_threshold( nb_e, bt2bt3_t, bt2bt3rel_t, bt1bt2_t, bt1bt2rel_t ) )
			{
				q.push( simple_pair( 0, nb_e, vi ) );
			}
		}
	}
	//cout << "after init, q size: " << q.size() << "" << endl;
	//cout << "prune: preparation done. " << endl;

	//vts_to_debug.insert(424);
	for (unsigned vi = 0; vi < m_vts.size(); ++vi) //debug
	{
		if (vts_to_debug.count(vi))
		{
			cout << "nbs of "<<vi<<": ";
			print_int_list(m_nb_edges_for_vert[vi]);
			cout << endl;
		}
	}

	prune_while_iteration(vts_to_debug, 
		bt2bt3_t, bt2bt3rel_t, bt1bt2_t, bt1bt2rel_t, q);

	// removed and ref count should be consistent
	/*cout << "sanity check - after pruning, before non-dual edge removal ..." << endl;
	sanity_check();*/

	/* remove non-dual edges (with ref-cnt == 0) & modify ref-cnt for vertices accordingly
	** Doing this might break topology, but removes unwanted connection,
	** e.g. the singular touching of 2 tri faces will be removed in the pruned skeleton
	** (appeared once in the neptune shape's MA)
	** However, it will still be preserved when using the new dual scheme. */
	for (unsigned ei = 0; ei < m_edges.size(); ++ei)
	{
		if (m_ref_edge_per_prune[ei] == 0 && part_of_orig_edge(ei) && !m_removed[0][ei])
		{
			m_removed[0][ei] = true;
			TriEdge e = m_edges[ei];
			m_ref_vert_per_prune[e[0]] --;
			m_ref_vert_per_prune[e[1]] --;

			// remove vert whose ref-cnt == 0
			// in case of both verts all having 0 ref-cnt, since preserving homotopy, 
			// at least one vert is remained
			int cnt = 0;
			if (m_ref_vert_per_prune[e[0]] == 0)
			{
				m_removed[2][e[0]] = true;
				cnt++;
			}
			if (m_ref_vert_per_prune[e[1]] == 0)
			{
				m_removed[2][e[1]] = true;
				cnt++;
			}
			if (cnt == 2)
				m_removed[2][e[1]] = false;				

			// debug
			if ( !(m_ref_vert_per_prune[e[0]] >= 0 && m_ref_vert_per_prune[e[1]] >= 0) )
			{
				cout << "logic err: ref-vert < 0!" << endl;
				exit(-1);
			}

			//cout << "non-dual edge removed." << endl;

			for (int i = 0; i < 2; ++i)
			{
				unsigned vi = e[i];
				if (m_ref_vert_per_prune[vi] == 1)
				{
					const auto& nbEdges = m_nb_edges_for_vert[vi];
					for (auto iter = nbEdges.begin(); iter != nbEdges.end(); ++iter)
					{
						if (!m_removed[0][*iter] && 
							edge_below_threshold(*iter, bt2bt3_t, bt2bt3rel_t, bt1bt2_t, bt1bt2rel_t))
						{
							q.push(simple_pair(0, *iter, vi));
							break;
						}
					}
				}
			}
		}
	}

	// removed and ref count should be consistent
	/*cout << "sanity check - after non-dual edge removal, before small components removal..." << endl;
	sanity_check();*/

	/* we may want to remove those small isolated components for better look */
	if (_remove_small_components)
	{
		cout << "removing *small* components..." << endl;

		// identify *small components* and remove all poly faces 
		// and non-dual edges, and all relevant vertices
		mark_components(
			m_removed, 
			m_ref_edge_per_prune, m_ref_vert_per_prune, 
			m_nb_faces_for_edge, m_nb_edges_for_vert, 
			m_face_degen_thresh, m_cmpnt_num_faces_thresh);
		//cout << "small components marked!" << endl;

		// need to re-initialize q with new simple pairs
		for (unsigned ei = 0; ei < m_edges.size(); ++ei)
		{
			if ( m_ref_edge_per_prune[ ei ] == 1 && !m_removed[ 0 ][ ei ] )
			{
				unsigned nb_f = m_nb_faces_for_edge[ei][0];
				const auto& nbFaces = m_nb_faces_for_edge[ei];
				for (auto nb_iter = nbFaces.begin(); nb_iter != nbFaces.end(); nb_iter++)
				{
					if ( !m_removed[1][*nb_iter] )
					{
						nb_f = *nb_iter;
						break;
					}
				}
				/*float val_bt2bt3 = edges_diff[ei];
				float val_bt2bt3rel = edges_reldiff[ei];*/
				if (face_edge_below_threshold(nb_f, ei, bt2bt3_t, bt2bt3rel_t, bt1bt2_t, bt1bt2rel_t))
				{
					q.push(simple_pair(1, nb_f, ei));
				}
			}
		}
	}

	// we may want to perform another pruning
	cout << "after marking small components, q size: " << q.size() << endl;
	if ( !q.empty() ) 
	{
		prune_while_iteration(vts_to_debug, 
			bt2bt3_t, bt2bt3rel_t, bt1bt2_t, bt1bt2rel_t, q);
	}

	// removed and ref count should be consistent
	/*cout << "sanity check - after small components removal..." << endl;
	sanity_check();*/

	/* remove non-dual edges (since we don't need them now anymore!) */
	//for (unsigned ei = 0; ei < m_edges.size(); ++ei)
	//{
	//	if (/*ref_edge[ei] == 0*/part_of_orig_edge(ei))
	//	{
	//		m_removed[0][ei] = true;
	//	}
	//}

	/* post pruning sanity check */
	// removed and ref count should be consistent
	/*cout << "sanity check - after (clean-up) non-dual edge removal..." << endl;
	sanity_check();*/

	// report what's left: any simple pair that should be removed however left?
	// first, edge-vert pair
	for (unsigned ei = 0; ei < m_edges.size(); ++ei)
	{
		if (m_removed[0][ei])
			continue;

		const auto& e = m_edges[ei];
		if ( m_ref_edge_per_prune[ei] == 1 && 
			edge_below_threshold(ei, bt2bt3_t, bt2bt3rel_t, bt1bt2_t, bt1bt2rel_t) )
		{
			if (m_ref_vert_per_prune[e[0]] == 1)
			{
				cout << "edge-vert pair "<< ei<<"-"<<e[0]<<"should be removed!"<<endl;
			}
			if (m_ref_vert_per_prune[e[1]] == 1)
			{
				cout << "edge-vert pair "<< ei<<"-"<<e[1]<<"should be removed!"<<endl;
			}
		}
	}
	// second, examine face-edge pair
	for (unsigned fi = 0; fi < m_tri_faces.size(); ++fi)
	{
		if (m_removed[1][fi])
			continue;

		const auto& f = m_tri_faces[fi];
		for (auto i = 0; i < 3; ++i)
		{
			unsigned ei = f[ i ];
			if ( m_ref_edge_per_prune[ei] == 1 && 
				!m_removed[0][ei] && !m_removed[1][fi] &&
				face_edge_below_threshold(fi, ei, bt2bt3_t, bt2bt3rel_t, bt1bt2_t, bt1bt2rel_t) )
			{
				cout << "face-edge pair "<<fi<<"-"<<ei<<"should have been removed!" <<endl
					<< "(check simple pair removal logic, or the logic used to perform this test.)"<<endl;
				//goto END_OF_POST_TEST;
			}
		}
	}
	// third, those faces marked as TO-REMOVE must have been removed
	for (unsigned fi = 0; fi < m_tri_faces.size(); ++fi)
	{
		if (m_to_remove_face[fi] && !m_removed[1][fi])
		{
			cout << "face " << fi << " (marked as to-remove) but not removed!"
				<< "value: " << face_bt2bt3[ fi ] << endl;
			const auto& f = m_tri_faces[fi];
			for ( auto i = 0; i < 3; ++i )
			{
				auto ei = f[ i ];
				if ( !m_removed[ 0 ][ ei ] )
				{
					cout << "edge " << ei << " removed? " << m_removed[ 0 ][ ei ] << ", "
						<< " dual edge? " << !part_of_orig_edge( ei ) << ", "
						<< " ref count " << m_ref_edge_per_prune[ ei ] << ", "
						<< " remaining nb faces: ";
					const auto& nbFaces = m_nb_faces_for_edge[ ei ];
					for ( auto f_iter = nbFaces.begin(); f_iter != nbFaces.end(); ++f_iter )
					{
						if ( !m_removed[ 1 ][ *f_iter ] )
							cout << *f_iter << ( m_to_remove_face[ *f_iter ] ? "(m)" : "" ) << ", ";
					}
					cout << endl;
				}
			}
			//goto END_OF_POST_TEST;
			//break;
		}
	}
END_OF_POST_TEST:
	;
}

void HybridSkeleton::exportSkeleton(
	std::string _skel_name, const trimesh::xform& _transform)
{
	// obtain remaining edges and faces of HS
	vector<TriEdge> edges_hs;
	getRemainedEdges(edges_hs);
	vector<TriFace> tri_faces_hs;
	getRemainedFaces(tri_faces_hs);
	vector<int> vts_remained_indices;
	getRemainedVts(vts_remained_indices);

	//cout << "# remaining vts ( >= # connected components): "<<vts_remained_indices.size()<<endl;

	// compact vertex set
	map<int, int> old_new_v_id_map;
	int new_id_cnt = 0;
	for (auto it = vts_remained_indices.begin(); it != vts_remained_indices.end(); ++it)
	{
		old_new_v_id_map[*it] = new_id_cnt;
		new_id_cnt ++;
	}
	// rename edge's vert ids
	for (auto it = edges_hs.begin(); it != edges_hs.end(); ++it)
	{
		auto& e = *it;
		e[0] = old_new_v_id_map.find(e[0])->second;
		e[1] = old_new_v_id_map.find(e[1])->second;
	}
	// rename face's vert ids
	for (auto it = tri_faces_hs.begin(); it != tri_faces_hs.end(); ++it)
	{
		auto& f = *it;
		f[0] = old_new_v_id_map.find(f[0])->second;
		f[1] = old_new_v_id_map.find(f[1])->second;
		f[2] = old_new_v_id_map.find(f[2])->second;
	}

	if ( _skel_name.find( ".sk" ) != std::string::npos )
	{
		// now export everything to file
		ofstream out_file( _skel_name );
		// # vts # edges # faces
		out_file << vts_remained_indices.size() << " " << edges_hs.size() << " " << tri_faces_hs.size() << endl;
		// vts coords...
		for ( auto it = vts_remained_indices.begin(); it != vts_remained_indices.end(); ++it )
		{
			auto v = m_vts[ *it ];
			v = _transform * v;
			out_file << v[ 0 ] << " " << v[ 1 ] << " " << v[ 2 ] << endl;
		}
		// edges...
		for ( auto it = edges_hs.begin(); it != edges_hs.end(); ++it )
		{
			const auto& e = *it;
			out_file << e[ 0 ] << " " << e[ 1 ] << endl;
		}
		// faces...
		for ( auto it = tri_faces_hs.begin(); it != tri_faces_hs.end(); ++it )
		{
			const auto& f = *it;
			out_file << f[ 0 ] << " " << f[ 1 ] << " " << f[ 2 ] << endl;
		}
		out_file.close();
	}
	else if ( _skel_name.find( ".ply" ) != std::string::npos )
	{
		// prepare vertices, edges, and faces for output
		vector<ply::Vertex> out_vts;
		vector<ply::Edge> out_edges;
		vector<ply::Face> out_faces;
		for ( auto it = vts_remained_indices.begin(); it != vts_remained_indices.end(); ++it )
		{
			auto v = m_vts[ *it ];
			v = _transform * v;
			out_vts.push_back( {v[0], v[1], v[2]} );
		}
		for ( auto it = edges_hs.begin(); it != edges_hs.end(); ++it )
		{
			const auto& e = *it;
			out_edges.push_back( { e[ 0 ], e[ 1 ] } );
		}
		for ( auto it = tri_faces_hs.begin(); it != tri_faces_hs.end(); ++it )
		{
			const auto& f = *it;
			out_faces.push_back( { 3, {f[ 0 ], f[ 1 ], f[ 2 ]} } );
		}
		std::map<std::string, PlyProperty> vert_props;
		vert_props[ "x" ] = { "x", Float32, Float32, offsetof( ply::Vertex, x ), PLY_SCALAR, 0, 0, 0 };
		vert_props[ "y" ] = { "y", Float32, Float32, offsetof( ply::Vertex, y ), PLY_SCALAR, 0, 0, 0 };
		vert_props[ "z" ] = { "z", Float32, Float32, offsetof( ply::Vertex, z ), PLY_SCALAR, 0, 0, 0 };
		std::map<std::string, PlyProperty> edge_props;
		edge_props[ "vertex1" ] = { "vertex1", Int32, Int32, offsetof( ply::Edge, v1 ), PLY_SCALAR, 0, 0, 0 };
		edge_props[ "vertex2" ] = { "vertex2", Int32, Int32, offsetof( ply::Edge, v2 ), PLY_SCALAR, 0, 0, 0 };
		std::map<std::string, PlyProperty> face_props;
		face_props[ "vertex_indices" ] = {
			"vertex_indices", Int32, Int32, offsetof( ply::Face, verts ),
			PLY_LIST, Uint8, Uint8, offsetof( ply::Face,nvts ) };
		ply::PLYwriter ply_writer;
		ply_writer.write( _skel_name.c_str(), true, true, true,
			vert_props, edge_props, face_props,
			out_vts, out_edges, out_faces );
		/*ply::PLYwriter ply_writer;
		ply_writer.write( _skel_name.c_str(), out_vts, out_edges, out_faces );*/
	}
}

void HybridSkeleton::absToRatio(
	float _abs_f_diff_t,
	float* _f_diff_r,
	float   _abs_f_reldiff_t,
	float* _f_reldiff_r,
	float   _abs_l_diff_t,
	float* _l_diff_r,
	float   _abs_l_reldiff_t,
	float* _l_reldiff_r
	)
{
	if (_f_diff_r)
		*_f_diff_r = (_abs_f_diff_t - min_bt2bt3) / (max_bt2bt3 - min_bt2bt3);
	if (_f_reldiff_r)
		*_f_reldiff_r = (_abs_f_reldiff_t - min_bt2bt3rel) / (max_bt2bt3rel - min_bt2bt3rel);
	if (_l_diff_r)
		*_l_diff_r = (_abs_l_diff_t - min_bt1bt2) / (max_bt1bt2 - min_bt1bt2);
	if (_l_reldiff_r)
		*_l_reldiff_r = (_abs_l_reldiff_t - min_bt1bt2rel) / (max_bt1bt2rel - min_bt1bt2rel);
}

void HybridSkeleton::ratioToAbs(
	float _f_diff_r /* = 0.0f */, float* _abs_f_diff_t /* = nullptr */, 
	float _f_reldiff_r /* = 0.0f */, float* _abs_f_reldiff_t /* = nullptr */, 
	float _l_diff_r /* = 0.0f */, float* _abs_l_diff_t /* = nullptr */, 
	float _l_reldiff_r /* = 0.0f */, float* _abs_l_reldiff_t /* = nullptr */)
{
	if (_abs_f_diff_t)
		*_abs_f_diff_t = min_bt2bt3 + (max_bt2bt3 - min_bt2bt3) * _f_diff_r;
	if (_abs_f_reldiff_t)
		*_abs_f_reldiff_t = min_bt2bt3rel + (max_bt2bt3rel - min_bt2bt3rel) * _f_reldiff_r;
	if (_abs_l_diff_t)
		*_abs_l_diff_t = min_bt1bt2 + (max_bt1bt2 - min_bt1bt2) * _l_diff_r;
	if (_abs_l_reldiff_t)
		*_abs_l_reldiff_t = min_bt1bt2rel + (max_bt1bt2rel - min_bt1bt2rel) * _l_reldiff_r;
}

void HybridSkeleton::to_vts_list( const TriFace & _tri_e_rep, int * _vts_l ) const
{
	auto e = m_edges[ _tri_e_rep[ 0 ] ];
	_vts_l[ 0 ] = e[ 0 ];
	_vts_l[ 1 ] = e[ 1 ];
	e = m_edges[ _tri_e_rep[ 1 ] ];
	_vts_l[ 2 ] = ( e[ 0 ] != _vts_l[ 1 ] ? e[ 0 ] : e[ 1 ] );
}

void HybridSkeleton::to_vts_list(const vector<unsigned>& _edge_l, vector<unsigned>& _vts_l) const
{
	_vts_l.reserve(_edge_l.size());
	auto start_e = m_edges[_edge_l[0]];
	auto second_e = m_edges[_edge_l[1]];
	// find the vertex to start from
	unsigned prev = start_e[0];
	if (second_e[0] == start_e[0] || second_e[1] == start_e[0])
		prev = start_e[1];
	_vts_l.push_back(prev);
	for (unsigned i = 0; i < _edge_l.size() - 1; ++i)
	{
		auto e = m_edges[_edge_l[i]];
		if (e[0] == prev)
		{
			_vts_l.push_back(e[1]);
			prev = e[1];
		}
		else
		{
			_vts_l.push_back(e[0]);
			prev = e[0];
		}
	}

	assert(_vts_l.size() == _edge_l.size());
}

bool HybridSkeleton::part_of_orig_edge(unsigned _ei) const
{
	return !m_stg->isDualEdgeInFineTri( transform_skel_edge( m_edges[ _ei ] ) );
}

int HybridSkeleton::transform_fine_vert_id( int _vi ) const
{
	if ( m_stg->getBurnScheme() == SteinerGraph::STEINER_ONLY )
	{
		bool is_dual = m_stg->isDualVertInFineTri( _vi, _vi );
		_vi += is_dual ? m_stg->m_stSubdiv.sizeOfStVts() : -m_stg->m_origG->vts.size();
	}
	return _vi;
}

TriEdge HybridSkeleton::transform_isolate_edge( const TriEdge & _e ) const
{
	if ( m_stg->getBurnScheme() == SteinerGraph::STEINER_ONLY )
		return _e + TriEdge( (int)m_stg->m_stSubdiv.sizeOfStVts() );
	return _e;
}

int HybridSkeleton::transform_skel_vert_id( int _vi ) const
{
	if ( m_stg->getBurnScheme() == SteinerGraph::STEINER_ONLY )
	{
		auto n_orig_vts = m_stg->m_origG->vts.size();
		_vi += n_orig_vts;
	}
	return _vi;
}

bool HybridSkeleton::edge_below_threshold(
	unsigned _ei, 
	float _bt2bt3_t, float _bt2bt3rel_t, float _bt1bt2_t, float _bt1bt2rel_t)
{
	if ( part_of_orig_edge(_ei) )
	{
		return true;
		/*if (edges_diff[_ei] < _bt2bt3_t || 
		edges_reldiff[_ei] < _bt2bt3rel_t)
		{
		return true;
		}*/
	}
	else
	{
		if ( edges_measure[ _ei ] < _bt1bt2_t ||
			edges_rel_measure[_ei] < _bt1bt2rel_t )
		{
			return true;
		}
	}
	return false;
}

bool HybridSkeleton::face_edge_below_threshold(
	unsigned _fi, unsigned _ei, 
	float _bt2bt3_t, float _bt2bt3rel_t, float _bt1bt2_t, float _bt1bt2rel_t)
{
	// if the face marked as removed & the edge is a non-dual edge, 
	// then remove the simple pair
	if (m_to_remove_face[_fi] && part_of_orig_edge(_ei))
	{
		//cout << "face-edge pair: face marked as removed; orig edge;" << endl;
		return true;
	}

	// else, only continue testing if it's below thresholds
	// always remove face-origedge pair
	// never remove face-dualedge pair
	if ( face_bt2bt3[_fi] < _bt2bt3_t/* ||
		face_bt2bt3rel[_fi] < _bt2bt3rel_t*/ )
	{
		return part_of_orig_edge( _ei );
	}
	return false;
}

void HybridSkeleton::prune_while_iteration(
	const set<unsigned>& vts_to_debug,
	float bt2bt3_t, float bt2bt3rel_t,
	float bt1bt2_t, float bt1bt2rel_t,
	queue<simple_pair>& q
	)
{
	int ei_debug = -1; // debug
	while ( !q.empty() )
	{
		const auto spr = q.front();
		q.pop();

		//debug
		if ( spr.type == 0 && vts_to_debug.count( spr.idx1 ) )
		{
			cout << "simple pair popped: (" << spr.type << ", " << spr.idx0 << ", " << spr.idx1 << ")"
				<< "ref_vert: " << m_ref_vert_per_prune[ spr.idx1 ] << endl;
		}
		if ( spr.type == 1 && spr.idx1 == ei_debug )
		{
			cout << "got face-edge spr (" << spr.idx0 << "," << spr.idx1 << ")"
				<< " face removed?" << ( m_removed[ 1 ][ spr.idx0 ] ? "true" : "false" ) << endl;
			print_edge( ei_debug );
		}

		// skip if this simple pair is removed already (the simple cell removed)
		if ( m_removed[ spr.type ][ spr.idx0 ] )
		{
			assert( spr.type == 1 ? m_ref_edge_per_prune[ spr.idx1 ] == 0 : m_ref_vert_per_prune[ spr.idx1 ] == 0 );

			//debug
			if ( spr.type == 0 && vts_to_debug.count( spr.idx1 ) )
			{
				cout << "skip: simple cell removed." << endl;
			}

			continue;
		}

		// logic checkpoint
		if ( !( spr.type == 1 && m_ref_edge_per_prune[ spr.idx1 ] == 1 || spr.type == 0 && m_ref_vert_per_prune[ spr.idx1 ] == 1 ) )
		{
			cout << "logic err: simple pair type " << spr.type;
			if ( spr.type == 1 )
				cout << " ref_edge[" << spr.idx1 << "] = " << m_ref_edge_per_prune[ spr.idx1 ];
			else
				cout << " ref_vert[" << spr.idx1 << "] = " << m_ref_vert_per_prune[ spr.idx1 ];
			cout << endl;
			exit( 1 );
		}

		// still valid simple pair. remove it, 
		// and modify ref count for the elements touched by the removal
		if ( spr.type == 0 ) // edge-vert pair
		{
			// set the simple pair: edge and vert removed
			m_removed[ 0 ][ spr.idx0 ] = true;
			m_removed[ 2 ][ spr.idx1 ] = true;

			// modify ref count for vts of the removed edge
			const auto& e = m_edges[ spr.idx0 ];

			m_ref_vert_per_prune[ e[ 0 ] ] --;
			m_ref_vert_per_prune[ e[ 1 ] ] --;

			//debug
			if ( vts_to_debug.count( spr.idx1 ) )
			{
				cout << "after dec (due to edge-vert removal "
					<< spr.idx0 << "-" << spr.idx1
					<< "), ref_vert["
					<< spr.idx1 << "]=" << m_ref_vert_per_prune[ spr.idx1 ] << endl;
			}

			assert( m_ref_vert_per_prune[ spr.idx1 ] == 0 );

			// add the only new simple pair edge-vert to q
			// if a new simple pair arises
			unsigned other_end = e[ 0 ] != spr.idx1 ? e[ 0 ] : e[ 1 ];
			if ( m_ref_vert_per_prune[ other_end ] == 1 )
			{
				const auto& nbs = m_nb_edges_for_vert[ other_end ];
				for ( auto nb_eit = nbs.begin(); nb_eit != nbs.end(); ++nb_eit )
				{
					if ( !m_removed[ 0 ][ *nb_eit ] )
					{
						if ( edge_below_threshold( *nb_eit, bt2bt3_t, bt2bt3rel_t, bt1bt2_t, bt1bt2rel_t ) )
						{
							q.push( simple_pair( 0, *nb_eit, other_end ) );

							if ( vts_to_debug.count( other_end ) ) //debug
							{
								cout << "simple pair added: ("
									<< 0 << ", " << *nb_eit << ", " << other_end << ")" << endl;
							}
						}
						break;
					}
				}
			}
		}
		else // face-edge pair
		{
			m_removed[ 1 ][ spr.idx0 ] = true;
			m_removed[ 0 ][ spr.idx1 ] = true;

			/*if (spr.idx0 == 176 && spr.idx1 == 1388)
			{
			cout << spr.idx0<<"removed? "<<removed[1][spr.idx0]<<" "
			<<spr.idx1<<"removed? "<<removed[0][spr.idx1]<<endl;
			}*/

			m_ref_edge_per_prune[ spr.idx1 ] --;
			assert( m_ref_edge_per_prune[ spr.idx1 ] == 0 );

			if ( spr.idx1 == ei_debug )
			{
				cout << "after face-edge spr (" << spr.idx0 << "," << spr.idx1 << ") removed" << endl;
				print_edge( ei_debug );
			}

			// modify ref count for the unremoved edges of the removed face
			const auto& f = m_tri_faces[ spr.idx0 ];
			for ( auto i = 0; i < 3; ++i )
			{
				auto ei = f[ i ];
				if ( !m_removed[ 0 ][ ei ] )
				{
					m_ref_edge_per_prune[ ei ] --;

					if ( spr.idx1 == ei_debug )
					{
						cout << "after modifying other edges of the face" << endl;
						print_edge( ei_debug );
					}

					// any new face-edge simple pair?
					if ( m_ref_edge_per_prune[ ei ] == 1 )
					{
						// find the unique face
						const auto& face_nbs = m_nb_faces_for_edge[ ei ];
						for ( auto nb_f = face_nbs.begin(); nb_f != face_nbs.end(); ++nb_f )
						{
							if ( !m_removed[ 1 ][ *nb_f ] )
							{
								// found the face. add to q the pair if criteria met
								if ( face_edge_below_threshold( *nb_f, ei, bt2bt3_t, bt2bt3rel_t, bt1bt2_t, bt1bt2rel_t ) )
								{
									q.push( simple_pair( 1, *nb_f, ei ) );
								}
								//if ( polyFace_bt2bt3[*nb_f] <= _bt2bt3_t ||
								//	polyFace_bt2bt3rel[*nb_f] <= _bt2bt3rel_t )
								//{
								//	if (part_of_orig_edge(*eit))
								//	{
								//		if (edges_diff[*eit] <= _bt2bt3_t || 
								//			edges_diff[*eit] <= _bt2bt3rel_t)
								//		{
								//			q.push(simple_pair(1, *nb_f, *eit));
								//		}
								//	}
								//	else
								//	{
								//		if (edges_diff[*eit] <= _bt1bt2_t || 
								//			edges_diff[*eit] <= _bt1bt2rel_t)
								//		{
								//			q.push(simple_pair(1, *nb_f, *eit));
								//		}
								//	}
								//}
								break;
							}
						}
					}
				}
			}

			// modify ref count for vts of the edge of the removed edge
			const auto& e = m_edges[ spr.idx1 ];
			for ( unsigned i = 0; i < 2; ++i )
			{
				unsigned vi = e[ (int)i ];
				m_ref_vert_per_prune[ vi ] --;
				if ( vts_to_debug.count( vi ) ) //debug
				{
					cout << "after dec (due to face-edge removal. type "
						<< spr.type << ", " << spr.idx0 << "-" << spr.idx1
						<< ", ref_vert[" << vi << "]=" << m_ref_vert_per_prune[ vi ] << endl;
				}

				// any new edge-vert simple pair?
				if ( m_ref_vert_per_prune[ vi ] == 1 )
				{
					// locate the only nb edge left for the vert
					const auto& nbs = m_nb_edges_for_vert[ vi ];
					for ( auto nb_ei = nbs.begin(); nb_ei != nbs.end(); ++nb_ei )
					{
						if ( !m_removed[ 0 ][ *nb_ei ] )
						{
							// edge located. add to q the pair if criteria met.
							if ( edge_below_threshold( *nb_ei, bt2bt3_t, bt2bt3rel_t, bt1bt2_t, bt1bt2rel_t ) )
							{
								q.push( simple_pair( 0, *nb_ei, vi ) );

								if ( vts_to_debug.count( vi ) ) //debug
								{
									cout << "simple pair added: ("
										<< 0 << ", " << *nb_ei << ", " << vi << ")" << endl;
								}
							}
							break;
						}
					}
				}
			}
		}
	}
	cout << "prune: iterative removal done. " << endl;
}

void HybridSkeleton::mark_components(
	const vector<bool>* const _removed, 
	const vector<int>& _ref_edge, const vector<int>& _ref_vert, 
	const vector<vector<unsigned>>& _nb_faces, const vector<vector<unsigned>>& _nb_edges,
	float _cmpnt_geomSize_t,
	int _num_face_t
	)
{
	vector<bool> edge_visited(m_edges.size(), false);

	/* 
	while there is unvisited edge:
		1. find a cmpnt C; mark edges as visited in the process;
		2. if |C| has too many faces, mark all its faces as NOT-TO-REMOVE;
		3. else, mark all faces as REMOVE, if C is *small* enough;
	*/
	queue<unsigned> q_faces; // queue of faces' id.s
	set<unsigned> face_cmpnt; // store faces of the cur component uniquely
	vector<TriPoint> vts_pos_cmpnt;
	vector<unsigned> edges_of_f;
	int cnt_total_cmpnts = 0;
	int cnt_small_cmpnts = 0;
	for (unsigned ei = 0; ei < m_edges.size(); ++ei)
	{
		if (_ref_edge[ei] != 1)
			continue;

		if (edge_visited[ei])
			continue;

		edge_visited[ei] = true;

		// start growing from this edge:
		// initialize queue with this only face
		const auto& nbFaces = _nb_faces[ei];
		int seed_fi = -1;
		for (unsigned i = 0; i < nbFaces.size(); ++i)
		{
			if ( !_removed[1][nbFaces[i]] )
			{
				seed_fi = nbFaces[i];
				break;
			}
		}
		assert(seed_fi >= 0);
		face_cmpnt.clear();
		face_cmpnt.insert(seed_fi);
		q_faces.push(seed_fi);

		// start growing from this face:
		edges_of_f.clear();
		while ( !q_faces.empty() )
		{
			auto fi = q_faces.front();
			q_faces.pop();

			const auto& f = m_tri_faces[fi];

			// find each nb face of this face f (a list of edge indices)
			for ( unsigned i = 0; i < 3/*f.size()*/; ++i )
			{
				if ( edge_visited[f[i]] )
					continue;

				edge_visited[f[i]] = true;
				// only inspect edges that are still used by other face(s)
				//if ( _ref_edge[f[i]] > 1 )
				{
					const auto& nbfaces_of_e = _nb_faces[f[i]];
					for (auto nb_iter = nbfaces_of_e.begin(); nb_iter != nbfaces_of_e.end(); ++nb_iter)
					{
						// skip nb face that's already removed
						if (_removed[1][*nb_iter] || *nb_iter == fi)
							continue;

						// only push those nb faces that are not yet in the current component to the queue 
						auto pre_size = face_cmpnt.size();
						face_cmpnt.insert(*nb_iter);
						q_faces.push(*nb_iter);
						//if ( pre_size < face_cmpnt.size() ) // face was not in the cmpnt?
						//{
						//	q_faces.push(*nb_iter);
						//}
					}
				}
			} // for (): done grabbing nb faces of f to component and queue
		} // while(): done for cur component

		cnt_total_cmpnts ++;

		// now component is ready. let's find out if it can be pruned.
		// is it too big?
		if (face_cmpnt.size() > _num_face_t)
		{
			// component is too big (too many faces, can't be an isolated small component, right?)
			// don't remove it. 
			continue;
		}
		else
		{
			// test if its span in the space is small enough (area)
			float sum_area = .0f;
			for (auto fi_iter = face_cmpnt.begin(); fi_iter != face_cmpnt.end(); ++fi_iter)
			{
				const auto& f = m_tri_faces[*fi_iter];
				sum_area += m_face_area[*fi_iter];
			}
			
			// is the component geometrically small enough?
			/*vts_pos_cmpnt.clear();
			for (auto iter = vts_cmpnt.begin(); iter != vts_cmpnt.end(); ++iter)
			{
			vts_pos_cmpnt.push_back(m_vts[*iter]);
			}
			bool is_geom_small = util::is_degenerate(vts_pos_cmpnt, _cmpnt_geomSize_t);*/
			
			bool is_geom_small = sum_area <= _cmpnt_geomSize_t;

			if (is_geom_small)
			{
				cnt_small_cmpnts ++;

				// this component is geometrically small. 
				// mark all its faces as TO-REMOVE!
				for (auto fi_iter = face_cmpnt.begin(); fi_iter != face_cmpnt.end(); ++fi_iter)
				{
					m_to_remove_face[*fi_iter] = true;
				}
			}
		} // if... else... : done testing whether to keep this component or not.
	} // for each edge: done looping thru each edge to locate disconnected components

	cout << cnt_total_cmpnts << " total components found. "<<endl;
	cout << cnt_small_cmpnts << " small components found. "<<endl;

} // mark_components()

void HybridSkeleton::print_edge( int _ei ) const
{
	const auto& e = m_edges[ _ei ];
	cout << "debug ---" << endl;
	cout << "edge: " << _ei << ", " << e << ". "
		<< "removed? "<<(m_removed[0][_ei] ? "true" : "false")<<". "
		<< "#ref-by-face " << m_ref_edge_per_prune[ _ei ] << ". " << endl;
	cout << "--- debug " << endl;
}

//
// debug & test routines
//
void print_direct_g(const map< unsigned, list<out_item> >& _directed_graph)
{
	for (auto it = _directed_graph.begin(); it != _directed_graph.end(); ++it)
	{
		cout << it->first <<": ";
		const auto& out_nodes = it->second;
		for (auto nit = out_nodes.begin(); nit != out_nodes.end(); ++nit)
		{
			cout << "("<<std::get<0>(*nit)<< "," <<std::get<1>(*nit)<<","<<std::get<2>(*nit)<<") -> ";
		}
		cout << '\b'<<'\b'<<'\b';
		cout << endl;
	}
}

void print_int_list(const vector<unsigned>& _l)
{
	for (auto it = _l.begin(); it != _l.end(); ++it)
		cout << *it<<" ";
}