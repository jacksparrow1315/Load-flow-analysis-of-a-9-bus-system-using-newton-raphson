clc;
clear;
close all;

nbus = input('Enter bus number-->');
Y = Y_bus(nbus);
busd = Bus_Data(nbus);
baseMVA = 100;

Newton_Raphson(nbus,Y,busd,baseMVA);