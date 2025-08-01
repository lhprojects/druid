#ifndef EVENT_NAV_
#define EVENT_NAV_
class EventNavigator
{
	public:
		void Fwd();
		void Bck();
		void setCollection();
		void MultiViewSwitch();
		void PTCutModify();
		void CellECutModify();
		void setCellColour(int ds);
		void HideLowSimE();
		void setPFOCellColour(int ds);
		void colorReroll();
		void SimuHitSizeModify();
		void PFOHitSizeModify();
		void HidePFOClu();
		void ClusterHitSizeModify();
		void ScaleModify();
		void HitTextAttach();
		void SizeModify();
		void setEnergyScale();
		void GotoEvent();
		void setGlobalEnergyScale();
		void setColorOverflowLimit();
		void setColorUnderflowLimit();
};
#endif //EVENT_NAV_

