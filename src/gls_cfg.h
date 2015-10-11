// gls_cfg.h: interface for the GLS_cfg class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(GLS_CFG_H__INCLUDED_)
#define GLS_CFG_H__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "fm_sys.h"

class GLS_cfg: public CFmSys
{
public:
    GLS_cfg();
    GLS_cfg(char* szCfgFile);
    virtual ~GLS_cfg();

    int m_nSectId;
    int m_nTaskId;
    int m_nMcmcHint;
    int m_nMcmcIter;
    int m_nMcmcSnps;
    bool m_bDebug;
    bool m_bRData;
    int m_lpOrder;
    int m_neffect;

    double m_fRhoTuning;
    double m_fDetRefit;
    double m_fPerRefit;
    double m_fPerSig;
    double m_fNumSig;
    double m_fQval_add;
    double m_fQval_dom;
    double m_fBurnInRound;

    int GetBurnInRound();

protected:
    int Load(char* szCfgFile);

};

void destroy(GLS_cfg* p);

#endif
