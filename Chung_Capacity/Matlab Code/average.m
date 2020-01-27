function average(x)
    figure();
    plot(x, x,'r');
    hold on
    plot(x, log2(x),'g');
    grid on;
end