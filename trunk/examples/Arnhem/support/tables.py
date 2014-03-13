x, y, corine_water, FillVal   = readMap(r'corine_water_sel_recl5.tif','GTiff')
corine_water[corine_water==FillVal]=0
x, y, corine_paved, FillVal   = readMap(r'corine_paved_sel_recl5.tif','GTiff')
corine_paved[corine_paved==FillVal]=0
x, y, GBKN_huizen, FillVal   = readMap(r'GBKN_gebouwen.tif','GTiff')
corine_water[corine_water==FillVal]=0
x, y, GBKN_wegen, FillVal   = readMap(r'GBKN_wegen.tif','GTiff')
corine_paved[corine_paved==FillVal]=0
x, y, GBKN_water, FillVal   = readMap(r'GBKN_water.tif','GTiff')
corine_paved[corine_paved==FillVal]=0
x, y, zones, FillValZones = readMap(r'buurten_project4.tif','GTiff')
icount=0    
output=np.zeros([10,len(list(set(zones.ravel())))])


for izone in list(set(zones.ravel())):
	output[0,icount]=izone
	print izone
	print icount
	print lu_roads_total.shape
	print output.shape
	output[1,icount]=np.mean(lu_roads_total[zones == izone])
	output[2,icount]=np.mean(lu_unpaved[zones == izone])
	output[3,icount]=np.mean(lu_buildings[zones == izone])
	output[4,icount]=np.mean(lu_water[zones == izone])
	output[5,icount]=np.mean(corine_water[zones == izone])
	output[6,icount]=np.mean(corine_paved[zones == izone])
	output[7,icount]=np.mean(GBKN_huizen[zones == izone])
	output[8,icount]=np.mean(GBKN_wegen[zones == izone])
	output[9,icount]=np.mean(GBKN_water[zones == izone])                
	icount=icount+1

np.savetxt("zones.csv", np.transpose(output), delimiter=",") 