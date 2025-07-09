import React from "react";
import { useNavigate } from "react-router-dom";
import { MdKeyboardDoubleArrowLeft } from "react-icons/md";

interface HomeButtonProps {
  className?: string;
  style?: React.CSSProperties;
  color?: string;
  hoverColor?: string;
}

const HomeButton: React.FC<HomeButtonProps> = ({ className = '', style = {}, color = '#16a34a', hoverColor }) => {
  const navigate = useNavigate();
  const [isHover, setIsHover] = React.useState(false);
  return (
    <button
      onClick={() => navigate("/")}
      className={`inline-flex items-center justify-center rounded-full focus:outline-none ${className}`}
      style={{ background: 'none', cursor: 'pointer', width: '44px', height: '44px', minWidth: '44px', minHeight: '44px', padding: 0, display: 'flex', alignSelf: 'flex-start', ...style }}
      title="Go to homepage"
      onMouseEnter={() => setIsHover(true)}
      onMouseLeave={() => setIsHover(false)}
    >
      <MdKeyboardDoubleArrowLeft size={40} color={hoverColor && isHover ? hoverColor : color} style={{ background: 'none', backgroundColor: 'transparent' }} />
    </button>
  );
};

export default HomeButton; 