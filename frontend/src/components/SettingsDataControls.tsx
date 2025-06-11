import React from "react";

const SettingsDataControls: React.FC = () => {
  return (
    <div className="flex items-center justify-center min-h-screen min-w-full w-full h-screen" style={{ background: 'linear-gradient(135deg, #1a1426 0%, #2a1a2a 60%, #231a23 100%)' }}>
      <div className="w-full max-w-2xl p-8 rounded shadow-lg" style={{ background: 'rgba(0,0,0,0.15)' }}>
        <h1 className="text-3xl font-bold mb-8 text-white text-center">Data Controls</h1>
        <div className="text-white text-center">This is the Data controls settings page.</div>
      </div>
    </div>
  );
};

export default SettingsDataControls; 