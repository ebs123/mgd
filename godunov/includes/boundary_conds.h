#pragma once

class CBoundaryConds
{
public:
	CBoundaryConds(void);
public:
	~CBoundaryConds(void);

private:
	void linear(double dl, double ul, double pl, double dr, double ur, double pr, double &d, double &u, double &p);
};
