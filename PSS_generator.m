function d_u= PSS_generator()  

u_shift = [25 29 34];
NID = 0;
d_u = [];

for n = 0:61
    u = u_shift(NID+1);
    if n <= 30
         d = exp(-j*pi*u*n*(n+1)/63);    
    else
         d = exp(-j*pi*u*(n+1)*(n+2)/63);   
    end;
    d_u = [d_u d];
end;

% figure
% subplot(1,3,1);
% plot(real(d_u(1:31)),imag(d_u(1:31)),'ko','MarkerFaceColor',[0 0 0]);
% axis([-1.5 1.5 -1.5 1.5]);
% title('n=0..30');
% subplot(1,3,2);
% plot(real(d_u(32:62)),imag(d_u(32:62)),'bo','MarkerFaceColor',[0 0 1]);
% axis([-1.5 1.5 -1.5 1.5]);
% title('n=31..61');
% subplot(1,3,3);
% plot(real(d_u(1:62)),imag(d_u(1:62)),'ro','MarkerFaceColor',[1 0 0]);
% axis([-1.5 1.5 -1.5 1.5]);
% title('n=0..61');
