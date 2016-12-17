pro readdata
DEVICE, DECOMPOSE = 0
loadct, 13

TVLCT, 255,255,255,254
TVLCT, 0,0,0,253
!P.Color = 253
!P.Background = 254

file1 = 'plot1.txt'
;file2 = 'plot2.txt'
;file3 = 'plot3.txt'
;file4 = 'plot4.txt'
;file5 = 'plot5.txt'

readcol, file1, tau,x,T,chi,p_tot
;readcol, file2, tau,P_gas,P_electron,P_ratio
;readcol, file3, tau,f11,f21,f12,f22,f32,f13,f23,f33
;readcol, file4, tau,n_electron,partial
;readcol, file5, tau,exc1,exc2,exc3

;x=alog(x) & tau=alog(tau) & T=alog(T) & chi=alog(chi) & p_tot=alog(p_tot)

;plot, tau,x,xtitle='optical depth',ytitle='log x(units)'
;layout=[2,2,1]
plot, tau,T,xtitle='optical depth',ytitle='log T(K)'
;layout=[2,2,2]
;plot, tau,chi,xtitle='optical depth',ytitle='Opacity',layout=[2,2,3]
;plot, tau,p_tot,xtitle='optical depth',ytitle='Total Mass Density ()',layout=[2,2,4]

end 
