; FOR WINDOWS ##### THIS MAY WORK ######
; LOAD RAI; SET UP HEATING
; BASICS col=2 is blue, col=4 is green, and col = 6 is red.

device, decomp=0
loadct, 39
!p.background = 255
!p.color = 0
tvlct, red, green, blue, /get

red(1) = 128
blue(1) = 128
green(1) = 0

red(2) = 0
blue(2) = 128
green(2) = 0

blue(3) = 128
red(3) = 0
green(3) = 128

blue(4) = 0
red(4) = 0
green(4) = 135

blue(5) = 0
red(5) = 200
green(5) = 175

blue(6) = 0
red(6) = 225
green(6) = 0

tvlct, red, green, blue
end
