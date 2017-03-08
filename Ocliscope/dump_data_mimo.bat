e:
cd E:\SoraSDK2.0\bin
dut start --radio 0
::Rx sample rate - 40MHz
dut radwr  --reg 0x17 --value 0x1 --radio 0

dut rxgain --value 0x800 --radio 0
dut rxpa --value 0x3000 --radio 0
dut centralfreq --value 5500 --radio 0
dut dump
echo File Dumped
pause

