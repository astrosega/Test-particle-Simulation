function finesfor, px, py, pz, fines

for i=0, n_elements(fines)-1 do begin
  
  if I eq 0 then a = px[fines[0],0]^2 + py[fines[0],0]^2 + pz[fines[0],0]^2 else a = [a, px[fines[j],j]^2 + py[fines[j],j]^2 + pz[fines[j],j]^2]
endfor

return, mean(a)
end