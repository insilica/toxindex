import React from "react";
import HomeButton from '../shared/HomeButton';

const SettingsDataControls: React.FC = () => {
  return (
    <div className="max-w-6xl mx-auto p-6 relative" style={{ paddingLeft: '8rem' }}>
      <div className="w-full max-w-2xl p-8 rounded shadow-lg" style={{ background: 'rgba(0,0,0,0.15)' }}>
        <h1 className="text-3xl font-bold mb-8 text-white text-center">Data Controls</h1>
        <div className="text-white text-center">This is the Data controls settings page.</div>
      </div>
      
      {/* Home button positioned relative to settings content */}
      <div 
        className="absolute transition-all duration-300 z-50"
        style={{
          left: '2rem',
          top: '1.5rem',
          border: 'none',
          padding: 0
        }}
      >
        <HomeButton
          color="#16a34a"
          hoverColor="#2563eb"
          aria-label="Go back"
        />
      </div>
    </div>
  );
};

export default SettingsDataControls; 