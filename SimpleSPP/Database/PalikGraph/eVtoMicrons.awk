BEGIN {
	h = 6.63e-34
	c = 3e8
	e = 1.6e-19

	unit = 1e6 #to convert m into µm
}
{
	print unit * h*c/($1*e), $2
}