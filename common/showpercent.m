function showpercent(i,n)

global spi;

curi= floor(i/n*100);

if i == 1 | isempty(spi)
    spi= curi;
    fprintf(1,'%2d%%',spi);
    return;
end

if spi < curi
    spi= curi;
    fprintf(1,'\b\b\b%2d%%',spi);
end

return
