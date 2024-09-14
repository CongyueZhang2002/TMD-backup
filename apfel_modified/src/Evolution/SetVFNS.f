************************************************************************
*
*     SetVFNS.f:
*
*     This subroutine sets the VFNS as a default.
*
************************************************************************
      subroutine SetVFNS
*
      implicit none
*
      include "../commons/Evs.h"
*
      Evs   = "VF"
      InEvs = "done"
*
      return
      end


************************************************************************
*
*     SetNPVFNS.f:
*
*     This subroutine sets the Non-Perturbative version of VFNS.
*
************************************************************************
      subroutine SetNPVFNS(np)
*
      implicit none
      logical np
*
      include "../commons/Evs.h"
*
      NPVF   = np
      InNpvf = "done"
*
      return
      end