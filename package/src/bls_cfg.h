// bls_cfg.h: interface for the BLS_cfg class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(BLS_CFG_H__INCLUDED_)
#define BLS_CFG_H__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "fm_sys.h"
#include "bls_R.h"

class BLS_cfg: public CFmSys
{
public:
    BLS_cfg();
    BLS_cfg(char* szCfgFile);
    virtual ~BLS_cfg();

    int Load(char* szCfgFile);

    int m_nSectId;
    int m_nTaskId;
    int m_nMcmcHint;
    int m_nMcmcIter;
    int m_nMcmcSnps;
    bool m_bDebug;
    bool m_bRData;
    int m_lpOrder;
    int m_neffect;
    bool m_bIncZ;
    bool m_bIncX;

    double m_fRhoTuning;
    double m_fDetRefit;
    double m_fPerRefit;
    double m_fPerSig;
    double m_fNumSig;
    double m_fQval_add;
    double m_fQval_dom;
    double m_fBurnInRound;

    int GetBurnInRound();

};

void destroy( BLS_cfg* p);

#endif
