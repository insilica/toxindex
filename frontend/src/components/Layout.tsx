import React, { useEffect, useState, useRef } from "react";
import { useNavigate } from "react-router-dom";

const Layout: React.FC<{ children: React.ReactNode }> = ({ children }) => {
  const navigate = useNavigate();
  const [auth, setAuth] = useState<null | boolean>(null);
  const [sidebarOpen, setSidebarOpen] = useState(true);
  const [sidebarClosing, setSidebarClosing] = useState(false);
  const [profileOpen, setProfileOpen] = useState(false);
  const profileRef = useRef<HTMLDivElement>(null);
  const [user, setUser] = useState<{ email?: string } | null>(null);

  useEffect(() => {
    fetch("/api/me", { credentials: "include", cache: "no-store" })
      .then(res => {
        if (res.status === 200) setAuth(true);
        else setAuth(false);
      })
      .catch(() => setAuth(false));
  }, []);

  useEffect(() => {
    fetch("/api/me", { credentials: "include", cache: "no-store" })
      .then(res => res.ok ? res.json() : null)
      .then(data => setUser(data));
  }, []);

  useEffect(() => {
    function handleClick(e: MouseEvent) {
      if (profileRef.current && !profileRef.current.contains(e.target as Node)) {
        setProfileOpen(false);
      }
    }
    if (profileOpen) document.addEventListener('mousedown', handleClick);
    return () => document.removeEventListener('mousedown', handleClick);
  }, [profileOpen]);

  const handleLogout = () => {
    fetch("/api/logout", {
      method: "GET",
      credentials: "include",
    }).then(() => {
      navigate("/login");
    });
  };

  return (
    <div className="flex w-screen min-h-screen overflow-x-hidden" style={{ fontFamily: 'Inter, Arial, sans-serif', position: 'relative' }}>
      {/* Sidebar */}
      <div style={{ position: 'relative' }}>
        <aside
          className={`bg-black shadow-sm flex flex-col justify-between min-h-screen relative transition-all duration-300 ease-in-out ${sidebarOpen ? 'w-64 opacity-100' : 'w-0 opacity-0 pointer-events-none'}`}
          style={{ overflow: 'hidden', position: 'relative', zIndex: 40, padding: 0, margin: 0, borderRight: 'none' }}
        >
          {/* Close button */}
          {sidebarOpen && (
            <button
              className={`absolute transition-all duration-300 w-9 h-9 flex items-center justify-center rounded-full bg-gray-800 text-gray-300 hover:bg-gray-700 hover:text-green-400 shadow`}
              style={{
                top: '2rem',
                zIndex: 50,
                minWidth: '36px',
                minHeight: '36px',
                padding: 0,
                right: sidebarOpen && !sidebarClosing ? '1.5rem' : 'auto',
                left: (!sidebarOpen || sidebarClosing) ? '0.5rem' : 'auto',
                opacity: sidebarOpen || sidebarClosing ? 1 : 0,
                pointerEvents: sidebarOpen || sidebarClosing ? 'auto' : 'none',
                transition: 'all 0.6s cubic-bezier(0.4,0,0.2,1)'
              }}
              onClick={() => {
                setSidebarClosing(true);
                setTimeout(() => {
                  setSidebarOpen(false);
                  setSidebarClosing(false);
                }, 300);
              }}
              title="Collapse sidebar"
            >
              <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5" strokeLinecap="round" strokeLinejoin="round">
                <polyline points="15 18 9 12 15 6" />
              </svg>
            </button>
          )}
          <div className="flex items-center mb-6" style={{ marginTop: '1.7rem', marginLeft: '1.2rem' }}>
            <button
              onClick={() => navigate('/')}
              className="inline-flex items-center justify-center rounded-full focus:outline-none"
              style={{ background: 'black', cursor: 'pointer', width: '44px', height: '44px', minWidth: '44px', minHeight: '44px', padding: 0, display: 'flex', alignSelf: 'flex-start' }}
              title="Go to homepage"
            >
              <svg width="44" height="44" viewBox="0 0 44 44" fill="none" xmlns="http://www.w3.org/2000/svg">
                <rect width="44" height="44" rx="10" fill="black"/>
                <text x="22" y="32" textAnchor="middle" fontSize="28" fontWeight="bold" fill="#16a34a" fontFamily="Arial, sans-serif">T</text>
              </svg>
            </button>
          </div>
          {/* Settings button at the bottom above logout */}
          <div className="flex flex-col items-center mt-auto mb-2">
            <a href="/settings" className="flex items-center hover:text-green-400 text-white px-4 py-2 rounded transition">
              <svg className="mr-2" width="22" height="22" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
                <circle cx="12" cy="12" r="3" />
                <path d="M19.4 15a1.65 1.65 0 0 0 .33 1.82l.06.06a2 2 0 0 1-2.83 2.83l-.06-.06a1.65 1.65 0 0 0-1.82-.33 1.65 1.65 0 0 0-1 1.51V21a2 2 0 0 1-4 0v-.09a1.65 1.65 0 0 0-1-1.51 1.65 1.65 0 0 0-1.82.33l-.06.06a2 2 0 0 1-2.83-2.83l.06-.06a1.65 1.65 0 0 0 .33-1.82 1.65 1.65 0 0 0-1.51-1H3a2 2 0 0 1 0-4h.09a1.65 1.65 0 0 0 1.51-1 1.65 1.65 0 0 0-.33-1.82l-.06-.06a2 2 0 0 1 2.83-2.83l.06.06a1.65 1.65 0 0 0 1.82.33h.09A1.65 1.65 0 0 0 9 3.09V3a2 2 0 0 1 4 0v.09a1.65 1.65 0 0 0 1 1.51h.09a1.65 1.65 0 0 0 1.82-.33l.06-.06a2 2 0 0 1 2.83 2.83l-.06.06a1.65 1.65 0 0 0-.33 1.82v.09a1.65 1.65 0 0 0 1.51 1H21a2 2 0 0 1 0 4h-.09a1.65 1.65 0 0 0-1.51 1z" />
              </svg>
              <span>Settings</span>
            </a>
          </div>
          {/* Logout button removed from sidebar. User profile/logout will be in dashboard. */}
        </aside>
        {/* Open sidebar button */}
        {!sidebarOpen && !sidebarClosing && (
          <button
            className="fixed left-2 z-50 w-9 h-9 flex items-center justify-center rounded-full bg-gray-800 text-gray-300 hover:bg-gray-700 hover:text-green-400 shadow transition"
            style={{ top: '2rem', minWidth: '36px', minHeight: '36px', padding: 0, position: 'fixed' }}
            onClick={() => setSidebarOpen(true)}
            title="Open sidebar"
          >
            <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5" strokeLinecap="round" strokeLinejoin="round">
              <polyline points="9 18 15 12 9 6" />
            </svg>
          </button>
        )}
      </div>
      {/* Model selection dropdown at top left of main app */}
      <div
        className="absolute top-0 w-full flex items-start justify-start z-30"
        style={{
          left: sidebarOpen || sidebarClosing ? '16rem' : 0,
          transition: 'left 0.3s cubic-bezier(0.4,0,0.2,1)'
        }}
      >
        <div
          className=""
          style={{
            marginTop: '2rem',
            marginLeft: sidebarOpen || sidebarClosing ? '2rem' : '3.5rem'
          }}
        >
          <div className="relative inline-block align-middle">
            <select
              className="text-white text-lg px-2 py-1 rounded-full border-none bg-transparent appearance-none transition hover:bg-gray-200 hover:bg-opacity-40 focus:bg-gray-200 focus:bg-opacity-40 focus:outline-none pr-8"
              style={{ minWidth: 120, boxShadow: 'none', background: 'none', cursor: 'pointer' }}
              defaultValue="toxindex-rap"
            >
              <option value="toxindex-rap">ToxIndex RAP</option>
              <option value="gpt-4">GPT-4</option>
              <option value="gpt-3.5">GPT-3.5</option>
              <option value="llama-2">Llama-2</option>
              <option value="custom">Custom Model</option>
            </select>
            <span style={{ pointerEvents: 'none', position: 'absolute', right: 6, top: '50%', transform: 'translateY(-50%)', display: 'flex', alignItems: 'center' }}>
              <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="#888" strokeWidth="2.2" strokeLinecap="round" strokeLinejoin="round">
                <path d="M7 10l5 5 5-5" />
              </svg>
            </span>
          </div>
        </div>
      </div>
      {/* User profile button at top right, aligned with others */}
      <div style={{ position: 'absolute', top: '2rem', right: '2.5rem', zIndex: 40 }} ref={profileRef}>
        <button
          onClick={() => setProfileOpen(v => !v)}
          className="w-11 h-11 flex items-center justify-center rounded-full bg-gray-900 bg-opacity-60 text-white hover:bg-gray-800 focus:outline-none shadow"
          style={{ padding: 0 }}
          title="User profile"
        >
          <svg width="28" height="28" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth="2" className="text-gray-300">
            <circle cx="12" cy="8" r="4" />
            <path d="M4 20c0-2.5 3.5-4 8-4s8 1.5 8 4" />
          </svg>
        </button>
        {profileOpen && (
          <div className="absolute right-0 mt-2 w-56 bg-white rounded-lg shadow-lg py-2 text-gray-900 border border-gray-200">
            <div className="px-4 py-2 border-b border-gray-100 text-sm">
              <div className="font-medium">{user?.email || "Logged in"}</div>
            </div>
            <button
              onClick={handleLogout}
              className="w-full text-left px-4 py-2 hover:bg-gray-100 text-red-600 text-sm"
            >
              Logout
            </button>
          </div>
        )}
      </div>
      {/* Main Content */}
      <main className="flex-1 flex items-center justify-center h-screen">
        {children}
      </main>
    </div>
  );
};

export default Layout; 