import React from "react";
import { useNavigate } from "react-router-dom";
import { MdKeyboardDoubleArrowLeft } from "react-icons/md";

interface HomeButtonProps {
  className?: string;
  style?: React.CSSProperties;
}

const HomeButton: React.FC<HomeButtonProps> = ({ className = '', style = {} }) => {
  const navigate = useNavigate();
  return (
    <button
      onClick={() => navigate("/")}
      className={`inline-flex items-center justify-center rounded-full focus:outline-none ${className}`}
      style={{ background: 'black', cursor: 'pointer', width: '44px', height: '44px', minWidth: '44px', minHeight: '44px', padding: 0, display: 'flex', alignSelf: 'flex-start', ...style }}
      title="Go to homepage"
    >
      {/* <svg width="44" height="44" viewBox="0 0 44 44" fill="none" xmlns="http://www.w3.org/2000/svg">
        <rect width="44" height="44" rx="10" fill="black"/>
        <text x="22" y="32" textAnchor="middle" fontSize="28" fontWeight="bold" fill="#16a34a" fontFamily="Arial, sans-serif">T</text>
      </svg> */}
      <MdKeyboardDoubleArrowLeft size={40} color="#16a34a" />
    </button>
  );
};

export default HomeButton; 