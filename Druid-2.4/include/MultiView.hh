#ifndef MULTIVIEW_
#define MULTIVIEW_

#include "TEveElement.h"
#include <TEveProjectionManager.h>
#include <TEveScene.h>
#include <TEveViewer.h>

//struct MultiView

class MultiView
{

	public:
	MultiView();

	~MultiView() {};

	void SetDepth(Float_t d);
	void ImportGeomRPhi(TEveElement* el);
	void ImportGeomRhoZ(TEveElement* el);
	void ImportEventRPhi(TEveElement* el);
	void ImportEventRhoZ(TEveElement* el);
	void DestroyEventRPhi();
	void DestroyEventRhoZ();

protected:
	TEveProjectionManager *fRPhiMgr;
	TEveProjectionManager *fRhoZMgr;
	TEveViewer            *f3DView;
	TEveViewer            *fRPhiView;
	TEveViewer            *fRhoZView;

	TEveScene             *fRPhiGeomScene;
	TEveScene             *fRhoZGeomScene;
	TEveScene             *fRPhiEventScene;
	TEveScene             *fRhoZEventScene;
};
#endif //MULTIVIEW_

