function f=fxy(x,y,fu)

if fu==1
    f=(((x.^2+y-11).^2+(x+y.^2-7).^2));
elseif fu==2
	f=20^2+x.^2+y.^2;
elseif fu==3
	f=x.^2./10^2-y.^2./15^2;
elseif fu==4
	f=cos(x).*y;
elseif fu==5
	f=cos(x).*cos(y).*exp(sqrt((x.^2+y.^2)./5));
end
