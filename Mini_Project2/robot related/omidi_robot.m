function [qnext]=omidi_robot(q,u,ts)
    qnext(1)=q(1)+u(1)*ts;
    qnext(2)=q(2)+u(2)*ts;

end
