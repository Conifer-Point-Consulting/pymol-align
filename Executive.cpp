#include"os_python.h"
#include"os_predef.h"
#include"os_std.h"
#include"os_gl.h"

#include <set>
#include <algorithm>
#include <map>
#include <cassert>
#include <clocale>
#include <vector>

#include "pymol/type_traits.h"

#ifdef PYMOL_OPENMP
#include <omp.h>
#endif

#include"Executive.h"
#include"Selector.h"
#include"Match.h"
#include"Util.h"
#include"Util2.h"
#include"PyMOL.h"
#include"PyMOLOptions.h"
#include"Parse.h"
#include"File.h"
#include"FileStream.h"

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#ifdef WIN32
#define mkstemp _mktemp_s
#endif

#define cTempRectSele "_rect"
#define cLeftButSele "lb"
#define cIndicateSele "indicate"

#ifndef NO_MMLIBS
#include "mmpymolx.h"
#endif


int ExecutiveAlign(PyMOLGlobals * G, const char *s1, const char *s2, const char *mat_file, float gap,
                   float extend, int max_gap, int max_skip, float cutoff, int cycles,
                   int quiet, const char *oname, int state1, int state2,
                   ExecutiveRMSInfo * rms_info, int transform, int reset, float seq_wt,
                   float radius, float scale, float base, float coord_wt, float expect,
                   int window, float ante)
{
  int sele1 = SelectorIndexByName(G, s1);
  int sele2 = SelectorIndexByName(G, s2);
  int *vla1 = NULL;
  int *vla2 = NULL;
  int na, nb;
  int c;
  int ok = true;
  int use_sequence = (mat_file && mat_file[0] && (seq_wt != 0.0F));
  int use_structure = (seq_wt >= 0.0F); /* negative seq_wt means sequence only! */
  ObjectMolecule *mobile_obj = NULL;
  CMatch *match = NULL;

  if(!use_structure)
    window = 0;

  if((scale == 0.0F) && (seq_wt == 0.0F) && (ante < 0.0F) && window)
    ante = window;

  if(ante < 0.0F)
    ante = 0.0F;

  if((sele1 >= 0)) {
    mobile_obj = SelectorGetSingleObjectMolecule(G, sele1);
    if(!mobile_obj) {
      ok = false;
      PRINTFB(G, FB_Executive, FB_Errors)
        " ExecutiveAlign: mobile selection must derive from one object only.\n" ENDFB(G);
    }
  }
  if(ok && (sele1 >= 0) && (sele2 >= 0) && rms_info) {
    vla1 = SelectorGetResidueVLA(G, sele1, use_structure, NULL);
    vla2 = SelectorGetResidueVLA(G, sele2, use_structure, mobile_obj);
    if(vla1 && vla2) {
      na = VLAGetSize(vla1) / 3;
      nb = VLAGetSize(vla2) / 3;
      if(na && nb) {
        match = MatchNew(G, na, nb, window);
        if(match) {
          if(use_sequence) {
            if(ok)
              ok = MatchResidueToCode(match, vla1, na);
            if(ok)
              ok = MatchResidueToCode(match, vla2, nb);
            if(ok)
              ok = MatchMatrixFromFile(match, mat_file, quiet);
            if(ok)
              ok = MatchPreScore(match, vla1, na, vla2, nb, quiet);
          }
          if(use_structure) {
	    /* avoid degenerate alignments */
	    ok = ((na>1) && (nb>1) && ok);
            if(ok) {
              ok = SelectorResidueVLAsTo3DMatchScores(G, match,
                                                      vla1, na, state1,
                                                      vla2, nb, state2, seq_wt,
                                                      radius, scale, base,
                                                      coord_wt, expect);
	    } else {
	      PRINTFB(G, FB_Executive, FB_Errors)
		" ExecutiveAlign: No alignment found.\n" ENDFB(G);
	    }
	  }
          if(ok)
            ok = MatchAlign(match, gap, extend, max_gap, max_skip, quiet, window, ante);
          if(ok) {
            rms_info->raw_alignment_score = match->score;
            rms_info->n_residues_aligned = match->n_pair;
            if(match->pair) {

              c = SelectorCreateAlignments(G, match->pair,
                                           sele1, vla1, sele2, vla2,
                                           "_align1", "_align2", false, false);

              if(c) {
                int mode = 2;
                if(!quiet) {
                  PRINTFB(G, FB_Executive, FB_Actions)
                    " %s: %d atoms aligned.\n", __func__, c ENDFB(G);
                }
                if(oname && oname[0] && reset)
                  ExecutiveDelete(G, oname);
                if(!transform)
                  mode = 1;
                ok = ExecutiveRMS(G, "_align1", "_align2", mode, cutoff, cycles,
                                  quiet, oname, state1, state2, false, 0, rms_info);
              } else {
                if(!quiet) {
                  PRINTFB(G, FB_Executive, FB_Actions)
                    " ExecutiveAlign-Error: atomic alignment failed (mismatched identifiers?).\n"
                    ENDFB(G);
                }
                ok = false;
              }
            }
          }
          MatchFree(match);
        }
      } else {
        ok = false;
        PRINTFB(G, FB_Executive, FB_Errors)
          " ExecutiveAlign: invalid selections for alignment.\n" ENDFB(G);
      }
    }
    ExecutiveUpdateCoordDepends(G, mobile_obj); //Updates dynamic_measures - see PYMOL-3090
  }

  VLAFreeP(vla1);
  VLAFreeP(vla2);
  return ok;
}


pymol::Result<std::vector<DiscardedRec>> ExecutiveDelete(PyMOLGlobals * G, pymol::zstring_view nameView, bool save)
{
  std::vector<DiscardedRec> discardedRecs;
  CExecutive *I = G->Executive;
  auto name = nameView.c_str();
  std::vector<OrderRec> specPositions;
  if (save) {
    specPositions = ExecutiveGetOrderOf(G, nameView);
  }

  auto getRecPos = [&] (SpecRec* rec) {
    auto it = std::find_if(specPositions.begin(), specPositions.end(),
                [&](OrderRec& oRec) {
                  return oRec.name == rec->name;
                });
    return it == specPositions.end() ? -1 : it->pos;
  };
  auto deleteObjRec = [&] (SpecRec* rec) {
    if(save) {
      if(rec->obj->type == cObjectGroup) {
        ExecutiveGroupPurge(G, rec, &discardedRecs);
      }
      ExecutivePurgeSpec(G, rec, save);
      auto rec_pos = getRecPos(rec);
      auto rec_out = ListDetachT(I->Spec, rec);
      SceneObjectRemove(G, rec->obj);
      assert(rec_pos);
      discardedRecs.emplace_back(rec_out, rec_pos);
    } else {
      if(rec->obj->type == cObjectGroup) {
        /* remove members of the group */
        ExecutiveGroup(G, rec->name, "", cExecutiveGroupPurge, true);
      }
      ExecutivePurgeSpec(G, rec, save);
      ListDelete(I->Spec, rec, next, SpecRec);        /* NOTE: order N in list length! TO FIX */
    }
  };

  auto deleteSeleRec = [&] (SpecRec* rec) {
    ExecutivePurgeSpec(G, rec, save);
    if (save) {
      auto rec_pos = getRecPos(rec);
      auto rec_out = ListDetachT(I->Spec, rec);
      assert(rec_pos);
      discardedRecs.emplace_back(rec_out, rec_pos);
    } else {
      ListDelete(I->Spec, rec, next, SpecRec);        /* NOTE: order N in list length! TO FIX */
    }
  };

  SpecRec *rec = NULL;
  CTracker *I_Tracker = I->Tracker;
  int list_id = ExecutiveGetNamesListFromPattern(G, name, false, false);
  int iter_id = TrackerNewIter(I_Tracker, 0, list_id);
  while(TrackerIterNextCandInList(I_Tracker, iter_id, (TrackerRef **) (void *) &rec)) {
    if(rec) {
      switch (rec->type) {
      case cExecSelection:
        deleteSeleRec(rec);
        break;
      case cExecObject:
        deleteObjRec(rec);
        break;
      case cExecAll:
        rec = nullptr;
        while(ListIterate(I->Spec, rec, next)) {
          switch (rec->type) {
          case cExecAll:
            break;
          case cExecObject:
            deleteObjRec(rec);
            rec = nullptr;
            break;
          case cExecSelection:
            deleteSeleRec(rec);
            rec = nullptr;
            break;
          }
        }
        SelectorDefragment(G);
        break;
      }
    }
  }
  TrackerDelList(I_Tracker, list_id);
  TrackerDelIter(I_Tracker, iter_id);

  // fix for PYMOL-757 - Note: This should probably go somewhere else, but we
  // couldn't figure out the correct place so far
  ExecutiveUpdateGroups(G, false);
  return discardedRecs;
}

int ExecutiveRMS(PyMOLGlobals * G, const char *s1, const char *s2, int mode, float refine,
                 int max_cyc, int quiet, const char *oname, int state1, int state2,
                 int ordered_selections, int matchmaker, ExecutiveRMSInfo * rms_info)
{
  /* mode 0: measure rms without superposition
     mode 1: measure rms with trial superposition
     mode 2: measure rms with actual superposition */

  int sele1, sele2;
  float rms = -1.0;
  int a, b;
  float inv, *f, *f1, *f2;
  ObjectMoleculeOpRec op1;
  ObjectMoleculeOpRec op2;
  OrthoLineType buffer;
  int *flag;
  int ok = true;
  int repeat;
  float v1[3], *v2;
  ObjectAlignment *align_to_update = NULL;

  bool ignore_case = SettingGetGlobal_b(G, cSetting_ignore_case);
  bool ignore_case_chain = SettingGetGlobal_b(G, cSetting_ignore_case_chain);

  int matrix_mode = SettingGetGlobal_i(G, cSetting_matrix_mode);
  if(matrix_mode < 0)
    matrix_mode = 0; /* for now */

  if(matchmaker == -1) {
    /* matchmaker -1 is the same as matchmaker 0 except that the
       selections are not pre-matched prior to calling of this routine */
    matchmaker = 0;
  }

  sele1 = SelectorIndexByName(G, s1);
  sele2 = SelectorIndexByName(G, s2);

  /* this function operates on stored coordinates -- thus transformation 
     matrices will need to be applied to the resulting atoms */

  // get coordinates
  {
    auto sele = sele1;
    auto state = state1;
    auto op = &op1;

    // for both selections
    do {
      ObjectMoleculeOpRecInit(op);

      if(sele >= 0) {
        if(state < 0) {
          op->code = OMOP_AVRT;
        } else {
          op->code = OMOP_StateVRT;
          op->i1 = state;
        }

        op->nvv1 = 0;                       // length of vc1 (number of atoms with coordinates)
        op->vc1 = VLACalloc(int, 1000);     // number of states per atom
        op->vv1 = VLACalloc(float, 1000);   // coordinates (sum over states)

        if(mode == 0)
          op->i2 = true;            /* if measuring current coordinates, then get global txfd values */

        if(matchmaker || (oname && oname[0]))
          op->ai1VLA = VLACalloc(AtomInfoType*, 1000);

        if(ordered_selections)
          op->vp1 = VLAlloc(int, 1000);     // selection member "priority"? (MemberType::tag)

        ExecutiveObjMolSeleOp(G, sele, op);

        for(a = 0; a < op->nvv1; a++) {
          inv = (float) op->vc1[a]; /* average over coordinate sets */
          if(inv) {
            f = op->vv1 + (a * 3);
            scale3f(f, 1.F / inv, f);
          }
        }
      }

      // second iteration done
      if (sele == sele2)
        break;

      sele = sele2;
      state = state2;
      op = &op2;

    } while (true);
  }

  if(op1.vv1 && op2.vv1) {
    if(op1.nvv1 && op2.nvv1) {
      ObjectMolecule *mobile_obj = NULL;

      int n_pair = 0;

      if(!(mobile_obj = SelectorGetSingleObjectMolecule(G, sele1))) {
        if(mode != 2) {
          PRINTFB(G, FB_Executive, FB_Warnings)
            "Executive-Warning: Mobile selection spans more than one object.\n" ENDFB(G);
        } else {
          PRINTFB(G, FB_Executive, FB_Errors)
            "Executive-Error: Mobile selection spans more than one object. Aborting.\n"
            ENDFB(G);
          ok = false;
        }
      }

      if(ok && op1.nvv1 && op2.nvv1 && (matchmaker > 0)) {
        /* matchmaker 0 is the default... by internal atom ordering only */
        int *idx1 = pymol::malloc<int>(op1.nvv1);
        int *idx2 = pymol::malloc<int>(op2.nvv1);
        int sort_flag = false;
        if(!(idx1 && idx2))
          ok = false;
        else {
          switch (matchmaker) {
          case 1:              /* by atom info-based ordering */
            UtilSortIndexGlobals(G, op1.nvv1, op1.ai1VLA, idx1,
                                 (UtilOrderFnGlobals *) fAtomOrdered);
            UtilSortIndexGlobals(G, op2.nvv1, op2.ai1VLA, idx2,
                                 (UtilOrderFnGlobals *) fAtomOrdered);
            sort_flag = true;
            break;
          case 2:              /* by matching atom identifiers */
            UtilSortIndex(op1.nvv1, op1.ai1VLA, idx1, (UtilOrderFn *) fAtomIDOrdered);
            UtilSortIndex(op2.nvv1, op2.ai1VLA, idx2, (UtilOrderFn *) fAtomIDOrdered);
            sort_flag = true;
            break;
          case 3:              /* by matching atom ranks */
            UtilSortIndex(op1.nvv1, op1.ai1VLA, idx1, (UtilOrderFn *) fAtomRankOrdered);
            UtilSortIndex(op2.nvv1, op2.ai1VLA, idx2, (UtilOrderFn *) fAtomRankOrdered);
            sort_flag = true;
            break;
          case 4:              /* by internal atom indexes (stored in temp1 kludge field) */
            UtilSortIndex(op1.nvv1, op1.ai1VLA, idx1, (UtilOrderFn *) fAtomTemp1Ordered);
            UtilSortIndex(op2.nvv1, op2.ai1VLA, idx2, (UtilOrderFn *) fAtomTemp1Ordered);
            sort_flag = true;
            break;
          }
          if(sort_flag) {
            /* GOD this is SO ugly! */

            if(op1.vv1) {
              float *tmp = VLAlloc(float, op1.nvv1 * 3);
              if(!tmp)
                ok = false;
              else {
                UtilApplySortedIndices(op1.nvv1, idx1, 3 * sizeof(float), op1.vv1, tmp);
                VLAFreeP(op1.vv1);
                op1.vv1 = tmp;
              }
            }
            if(op1.vc1) {
              int *tmp = VLAlloc(int, op1.nvv1);
              if(!tmp)
                ok = false;
              else {
                UtilApplySortedIndices(op1.nvv1, idx1, sizeof(int), op1.vc1, tmp);
                VLAFreeP(op1.vc1);
                op1.vc1 = tmp;
              }
            }
            if(op1.vp1) {
              int *tmp = VLAlloc(int, op1.nvv1);
              if(!tmp)
                ok = false;
              else {
                UtilApplySortedIndices(op1.nvv1, idx1, sizeof(int), op1.vp1, tmp);
                VLAFreeP(op1.vp1);
                op1.vp1 = tmp;
              }
            }
            if(op1.ai1VLA) {
              AtomInfoType **tmp = VLACalloc(AtomInfoType *, op1.nvv1);
              if(!tmp)
                ok = false;
              else {
                UtilApplySortedIndices(op1.nvv1, idx1, sizeof(AtomInfoType *), op1.ai1VLA,
                                       tmp);
                VLAFreeP(op1.ai1VLA);
                op1.ai1VLA = tmp;
              }
            }

            if(op2.vv1) {
              float *tmp = VLAlloc(float, op2.nvv1 * 3);
              if(!tmp)
                ok = false;
              else {
                UtilApplySortedIndices(op2.nvv1, idx2, 3 * sizeof(float), op2.vv1, tmp);
                VLAFreeP(op2.vv1);
                op2.vv1 = tmp;
              }
            }
            if(op2.vc1) {
              int *tmp = VLAlloc(int, op2.nvv1);
              if(!tmp)
                ok = false;
              else {
                UtilApplySortedIndices(op2.nvv1, idx2, sizeof(int), op2.vc1, tmp);
                VLAFreeP(op2.vc1);
                op2.vc1 = tmp;
              }
            }
            if(op2.vp1) {
              int *tmp = VLAlloc(int, op2.nvv1);
              if(!tmp)
                ok = false;
              else {
                UtilApplySortedIndices(op2.nvv1, idx2, sizeof(int), op2.vp1, tmp);
                VLAFreeP(op2.vp1);
                op2.vp1 = tmp;
              }
            }
            if(op2.ai1VLA) {
              AtomInfoType **tmp = VLACalloc(AtomInfoType *, op2.nvv1);
              if(!tmp)
                ok = false;
              else {
                UtilApplySortedIndices(op2.nvv1, idx2, sizeof(AtomInfoType *), op2.ai1VLA,
                                       tmp);
                VLAFreeP(op2.ai1VLA);
                op2.ai1VLA = tmp;
              }
            }

          }
        }

        if(matchmaker != 0) {
          int n1 = 0, n2 = 0, c1 = 0, c2 = 0;
          int cmp;

          while((n1 < op1.nvv1) && (n2 < op2.nvv1)) {
            cmp = 0;
            switch (matchmaker) {
            case 1:            /* insure that AtomInfoType matches */
              if(AtomInfoMatch(G, op1.ai1VLA[n1], op2.ai1VLA[n2], ignore_case, ignore_case_chain))
                cmp = 0;
              else
                cmp = AtomInfoCompare(G, op1.ai1VLA[n1], op2.ai1VLA[n2]);
              printf("%d-%d %d-%d: %d\n", c1, n1, c2, n2, cmp);
              break;
            case 2:            /* ID */
            case 3:            /* rank */
              {
                int val1;
                int val2;

                switch (matchmaker) {
                case 2:        /* ID */
                  val1 = op1.ai1VLA[n1]->id;
                  val2 = op2.ai1VLA[n2]->id;
                  break;
                case 3:        /* rank */
                  val1 = op1.ai1VLA[n1]->rank;
                  val2 = op2.ai1VLA[n2]->rank;
                  break;
                case 4:        /* index (via temp1) */
                  val1 = op1.ai1VLA[n1]->temp1;
                  val2 = op2.ai1VLA[n2]->temp1;
                  break;
                default:
                  val1 = 0;
                  val2 = 0;
                  break;
                }
                if(val1 == val2)
                  cmp = 0;
                else if(val1 < val2)
                  cmp = -1;
                else
                  cmp = 1;
              }
              break;
            }
            if(!cmp) {          /* match found */
              idx1[c1++] = n1++;
              idx2[c2++] = n2++;
              n_pair++;
            } else if(cmp < 0) {        /* op1 below op2 */
              n1++;
            } else {            /* op2 below op1 */
              n2++;
            }
          }

          if(n_pair) {
            if(op1.vv1)
              PackSortedIndices(n_pair, idx1, 3 * sizeof(float), op1.vv1);
            if(op1.vc1)
              PackSortedIndices(n_pair, idx1, sizeof(int), op1.vc1);
            if(op1.vp1)
              PackSortedIndices(n_pair, idx1, sizeof(int), op1.vp1);
            if(op1.ai1VLA)
              PackSortedIndices(n_pair, idx1, sizeof(AtomInfoType *), op1.ai1VLA);

            if(op2.vv1)
              PackSortedIndices(n_pair, idx2, 3 * sizeof(float), op2.vv1);
            if(op2.vc1)
              PackSortedIndices(n_pair, idx2, sizeof(int), op2.vc1);
            if(op2.vp1)
              PackSortedIndices(n_pair, idx2, sizeof(int), op2.vp1);
            if(op2.ai1VLA)
              PackSortedIndices(n_pair, idx2, sizeof(AtomInfoType *), op2.ai1VLA);
          }
        }
        FreeP(idx1);
        FreeP(idx2);
      } else if(op1.nvv1 != op2.nvv1) {
        sprintf(buffer, "Atom counts between selections don't match (%d vs %d)",
                op1.nvv1, op2.nvv1);
        ErrMessage(G, "ExecutiveRMS", buffer);
        n_pair = 0;
        ok = false;
      } else {
        n_pair = 0;
        for(a = 0; a < op1.nvv1; ++a) { // for atoms in selection
          if (op1.vc1[a] && op2.vc1[a]) { // check state counts
            if (n_pair < a) { // copy over if necessary
              copy3(op1.vv1 + 3 * a, op1.vv1 + 3 * n_pair);
              copy3(op2.vv1 + 3 * a, op2.vv1 + 3 * n_pair);
              if(op1.ai1VLA) op1.ai1VLA[n_pair] = op1.ai1VLA[a];
              if(op2.ai1VLA) op2.ai1VLA[n_pair] = op2.ai1VLA[a];
              if(op1.vp1) op1.vp1[n_pair] = op1.vp1[a];
              if(op2.vp1) op2.vp1[n_pair] = op2.vp1[a];
              op1.vc1[n_pair] = op1.vc1[a];
              op2.vc1[n_pair] = op2.vc1[a];
            }
            ++n_pair;
          }
        }
      }

      if(n_pair) {
        /* okay -- we're on track to do an alignment */

        if(ordered_selections && op1.vp1 && op2.vp1) {
          /* if we expected ordered selections and have priorities, 
             then we may need to sort vertices */

          int sort_flag1 = false, sort_flag2 = false;
          int well_defined1 = true, well_defined2 = true;

          for(a = 0; a < (n_pair - 1); a++) {
            /*          printf("op1 vertex %d priority %d\n",a,op1.vp1[a]);
               printf("op2 vertex %d priority %d\n",a,op2.vp1[a]); */

            if(op1.vp1[a] > op1.vp1[a + 1])
              sort_flag1 = true;
            else if(op1.vp1[a] == op1.vp1[a + 1])
              well_defined1 = false;
            if(op2.vp1[a] > op2.vp1[a + 1])
              sort_flag2 = true;
            else if(op2.vp1[a] == op2.vp1[a + 1])
              well_defined2 = false;
          }

          if(sort_flag1 || sort_flag2) {
            if(!(well_defined1 || well_defined2)) {
              PRINTFB(G, FB_Executive, FB_Warnings)
                "Executive-Warning: Ordering requested but not well defined.\n" ENDFB(G);
            } else {
              FitVertexRec *vert = pymol::malloc<FitVertexRec>(n_pair);

              if(sort_flag1) {
                float *src, *dst;
                src = op1.vv1;
                for(a = 0; a < n_pair; a++) {
                  vert[a].priority = op1.vp1[a];
                  dst = vert[a].vertex;
                  copy3f(src, dst);
                  src += 3;
                }
                UtilSortInPlace(G, vert, n_pair, sizeof(FitVertexRec),
                                (UtilOrderFn *) fVertexOrdered);
                dst = op1.vv1;
                for(a = 0; a < n_pair; a++) {
                  src = vert[a].vertex;
                  copy3f(src, dst);
                  dst += 3;
                }
              }

              if(sort_flag2) {
                float *src, *dst;
                src = op2.vv1;
                for(a = 0; a < n_pair; a++) {
                  vert[a].priority = op2.vp1[a];
                  dst = vert[a].vertex;
                  copy3f(src, dst);
                  src += 3;
                }
                UtilSortInPlace(G, vert, n_pair, sizeof(FitVertexRec),
                                (UtilOrderFn *) fVertexOrdered);
                dst = op2.vv1;
                for(a = 0; a < n_pair; a++) {
                  src = vert[a].vertex;
                  copy3f(src, dst);
                  dst += 3;
                }
              }

              FreeP(vert);
            }
          }
        }

        if(rms_info) {
          rms_info->initial_n_atom = n_pair;
          rms_info->n_cycles_run = 0;
          rms_info->final_n_atom = n_pair;      /* in case there is no refinement */
        }

        if(mode != 0) {
          rms = MatrixFitRMSTTTf(G, n_pair, op1.vv1, op2.vv1, NULL, op2.ttt);
          if(rms_info) {
            rms_info->initial_rms = rms;
            rms_info->final_rms = rms;
          }
          repeat = true;
          b = 0;
          while(repeat) {
            repeat = false;
            b++;
            if(b > max_cyc)
              break;
            if((refine > R_SMALL4) && (rms > R_SMALL4)) {
              int n_next = n_pair;
              AtomInfoType **ai1, **ai2;

              flag = pymol::malloc<int>(n_pair);

              if(flag) {
                for(a = 0; a < n_pair; a++) {
                  MatrixTransformTTTfN3f(1, v1, op2.ttt, op1.vv1 + (a * 3));
                  v2 = op2.vv1 + (a * 3);
                  if((diff3f(v1, v2) / rms) > refine) {
                    flag[a] = false;
                    repeat = true;
                  } else
                    flag[a] = true;
                }
                f1 = op1.vv1;
                f2 = op2.vv1;
                ai1 = op1.ai1VLA;
                ai2 = op2.ai1VLA;
                for(a = 0; a < n_pair; a++) {
                  if(!flag[a]) {
                    n_next--;
                  } else {
                    copy3f(op1.vv1 + (3 * a), f1);
                    copy3f(op2.vv1 + (3 * a), f2);
                    f1 += 3;
                    f2 += 3;
                    if(ai1 && ai2) {    /* make sure we keep track of which atoms are aligned */
                      *(ai1++) = op1.ai1VLA[a];
                      *(ai2++) = op2.ai1VLA[a];
                    }
                  }
                }
                if(!quiet && (n_next != n_pair)) {
                  PRINTFB(G, FB_Executive, FB_Actions)
                    " %s: %d atoms rejected during cycle %d (RMSD=%0.2f).\n", __func__,
                    n_pair - n_next, b, rms ENDFB(G);
                }
                n_pair = n_next;
                FreeP(flag);
                if(n_pair) {
                  rms = MatrixFitRMSTTTf(G, n_pair, op1.vv1, op2.vv1, NULL, op2.ttt);
                  if(rms_info) {
                    rms_info->n_cycles_run = b;
                    rms_info->final_n_atom = n_pair;
                    rms_info->final_rms = rms;
                  }
                } else
                  break;
              }
            }
          }
        } else {                /* mode == 0 -- simple RMS, with no coordinate movement */
          rms = MatrixGetRMS(G, n_pair, op1.vv1, op2.vv1, NULL);
          if(rms_info) {
            rms_info->initial_rms = rms;
            rms_info->final_rms = rms;
          }
        }
      }
      if(!n_pair) {
        PRINTFB(G, FB_Executive, FB_Results)
          " Executive: Error -- no atoms left after refinement!\n" ENDFB(G);
        ok = false;
      }

      if(ok) {
        if(!quiet) {
          PRINTFB(G, FB_Executive, FB_Results)
            " Executive: RMSD = %8.3f (%d to %d atoms)\n", rms, n_pair, n_pair ENDFB(G);
        }
        if(oname && oname[0]) {
            int align_state = state2;
            ObjectMolecule *trg_obj = SelectorGetSingleObjectMolecule(G, sele2);

            if(align_state < 0) {
              align_state = SceneGetState(G);
            }

            /* we're going to create/update an alignment object */

            {
              /* Get unique ids and construct the alignment vla */
              pymol::vla<int> align_vla(n_pair * 3);

              {
                int *id_p = align_vla.data();
                int i;
                for(i = 0; i < n_pair; i++) {
                  id_p[0] = AtomInfoCheckUniqueID(G, op2.ai1VLA[i]);    /* target */
                  id_p[1] = AtomInfoCheckUniqueID(G, op1.ai1VLA[i]);
                  id_p[2] = 0;
                  id_p += 3;
                }
                VLASize(align_vla, int, n_pair * 3);
              }
              {
                ObjectAlignment *obj = NULL;

                /* does object already exist? */
                {
                  pymol::CObject *execObj = ExecutiveFindObjectByName(G, oname);
                  if(execObj && (execObj->type != cObjectAlignment))
                    ExecutiveDelete(G, oname);
                  else
                    obj = (ObjectAlignment *) execObj;
                }
                obj =
                  ObjectAlignmentDefine(G, obj, align_vla, align_state, true, trg_obj,
                                        mobile_obj);
                obj->Color = ColorGetIndex(G, "yellow");
                ObjectSetName(obj, oname);
                ExecutiveManageObject(G, obj, 0, quiet);
                align_to_update = obj;
                SceneInvalidate(G);
              }
            }
        }
        if(ok && mode == 2) {
          if(matrix_mode>0) {

            ObjectMolecule *src_obj, *trg_obj;
            src_obj = SelectorGetFirstObjectMolecule(G, sele1); /* get at least one object */
            trg_obj = SelectorGetSingleObjectMolecule(G, sele2);

            /* first we need to make sure that the object being moved
               matches the target with respect to both the TTT and the
               object's state matrix (if any) */

            if(src_obj && trg_obj) {
              ExecutiveMatrixCopy(G, trg_obj->Name, src_obj->Name, 1, 1,        /* TTT mode */
                                  state2, state1, false, 0, quiet);

              ExecutiveMatrixCopy(G, trg_obj->Name, src_obj->Name, 2, 2,        /* Object state mode */
                                  state2, state1, false, 0, quiet);

              switch (matrix_mode) {
              case 1:          /* TTTs */
                ExecutiveCombineObjectTTT(G, src_obj->Name, op2.ttt, true, -1);
                break;
              case 2:
                {
                  double homo[16], *src_homo;
                  convertTTTfR44d(op2.ttt, homo);
                  if(ExecutiveGetObjectMatrix
                     (G, src_obj->Name, state1, &src_homo, false)) {
                    left_multiply44d44d(src_homo, homo);
                    ExecutiveSetObjectMatrix(G, src_obj->Name, state1, homo);
                  }
                }
                break;
              }
              /* next we need to update the object's TTT matrix to reflect
                 the transformation */
            }
          } else {              /* matrix_mode is zero -- legacy behavior */
            /* this will transform the actual coordinates */
            op2.code = OMOP_TTTF;
            ExecutiveObjMolSeleOp(G, sele1, &op2);
          }
        }
      }
    } else {
      ErrMessage(G, __func__, "No atoms selected.");
      ok = false;
    }
  }

  if(align_to_update) {
    align_to_update->update();
  }

  VLAFreeP(op1.vv1);
  VLAFreeP(op2.vv1);
  VLAFreeP(op1.vc1);
  VLAFreeP(op2.vc1);
  VLAFreeP(op1.vp1);
  VLAFreeP(op2.vp1);
  VLAFreeP(op1.ai1VLA);
  VLAFreeP(op2.ai1VLA);
  return (ok);
}

void ExecutiveUpdateCoordDepends(PyMOLGlobals * G, ObjectMolecule * mol)
{                               /* nasty, ugly, inefficient hack */

  CExecutive *I = G->Executive;
  SpecRec *rec = NULL;
  ObjectGadget *gadget;
  int done_inv_all = false;
  int dynamic_measures = SettingGet_b(G, mol ? mol->Setting.get() : NULL, NULL,
      cSetting_dynamic_measures);

  while(ListIterate(I->Spec, rec, next)) {
    if(rec->type == cExecObject) {
      switch(rec->obj->type) {
      case cObjectGadget:
        if(done_inv_all)
          break;
        gadget = (ObjectGadget *) rec->obj;
        if(gadget->GadgetType == cGadgetRamp) {
          ObjectGadgetRamp *ramp = (ObjectGadgetRamp *) gadget;
          if(ramp->RampType == cRampMol) {
            if(ramp->Mol == mol) {
              ExecutiveInvalidateRep(G, cKeywordAll, cRepAll, cRepInvColor);
              done_inv_all = true;
            }
          }
        }
        break;
      case cObjectMeasurement:
        if(dynamic_measures)
          ObjectDistMoveWithObject((ObjectDist*) rec->obj, mol);
        break;
      case cObjectAlignment:
        rec->obj->invalidate(
            cRepAll, cRepInvRep, cSelectorUpdateTableAllStates);
        break;
      }
    }
  }
}