# tell 
empirical min loop number
wehre/when might error be thrown

see where time checks are

currently plant-soil solute flow defined outside of iteration loop. will be different for N or P simulation

9c:
1		2	3	4	5	6	7	8		9
water,	cs,	cl,	b1,	b2,	b3,	b4,	css2,	co2

also called 10c because counts the soil even if it s static
by setting no microbes, can get normal richards2c simulation

check that actually get low SinkLim1DS, OutLim1DS, InOutBC_Cdiff

# to do
finish cleaning the cyl3plant and fpithelper.
test for a smal simulation of 1h
separate phloema and plant water flow, adapt transpiration
go through functions again to put unit and so on