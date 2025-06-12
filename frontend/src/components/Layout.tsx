import React, { useEffect, useState, useRef, isValidElement, cloneElement } from "react";
import { useNavigate, useLocation } from "react-router-dom";
import Dashboard from "./Dashboard";
import FilePreviewModal from './FilePreviewModal';
import { FaComments, FaPlus, FaListAlt } from 'react-icons/fa';

const Layout: React.FC<{ children: React.ReactNode }> = ({ children }) => {
  const navigate = useNavigate();
  const location = useLocation();
  const [auth, setAuth] = useState<null | boolean>(null);
  const [sidebarOpen, setSidebarOpen] = useState(true);
  const [sidebarClosing, setSidebarClosing] = useState(false);
  const [profileOpen, setProfileOpen] = useState(false);
  const profileRef = useRef<HTMLDivElement>(null);
  const [user, setUser] = useState<{ email?: string } | null>(null);
  const [selectedModel, setSelectedModel] = useState("toxindex-rap");
  const [environments, setEnvironments] = useState<{ environment_id: string; title: string }[]>([]);
  const [selectedEnv, setSelectedEnv] = useState<string>("");
  const [envFiles, setEnvFiles] = useState<{ file_id: number; filename: string }[]>([]);
  const [sidebarPreviewFileId, setSidebarPreviewFileId] = useState<number | null>(null);
  const [sidebarPreviewOpen, setSidebarPreviewOpen] = useState(false);
  const [chatSessions, setChatSessions] = useState<any[]>([]);
  const [selectedSessionId, setSelectedSessionId] = useState<string | null>(null);

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

  // Fetch files for selected environment
  useEffect(() => {
    if (!selectedEnv) return;
    fetch(`/api/environment/${selectedEnv}/files`, { credentials: 'include' })
      .then(res => res.json())
      .then(data => setEnvFiles(data.files || []));
  }, [selectedEnv]);

  // Fetch environments and set default selectedEnv
  useEffect(() => {
    fetch("/api/environments", { credentials: "include" })
      .then(res => res.json())
      .then(data => {
        setEnvironments(data.environments || []);
        if (data.environments && data.environments.length > 0) {
          setSelectedEnv(data.environments[0].environment_id);
        } else {
          setSelectedEnv("__add__");
        }
      });
  }, []);

  // Fetch chat sessions for environment
  useEffect(() => {
    if (!selectedEnv) return;
    fetch(`/api/environment/${selectedEnv}/chat_sessions`, { credentials: 'include' })
      .then(res => res.json())
      .then(data => {
        setChatSessions(data.sessions || []);
        if (data.sessions && data.sessions.length > 0) {
          setSelectedSessionId(data.sessions[0].session_id);
        } else {
          setSelectedSessionId(null);
        }
      });
  }, [selectedEnv]);

  const handleLogout = () => {
    fetch("/api/logout", {
      method: "GET",
      credentials: "include",
    }).then(() => {
      navigate("/login");
    });
  };

  // Settings sidebar logic
  const isSettings = location.pathname.startsWith("/settings");
  let settingsSection: 'general' | 'environments' | 'data-controls' = 'general';
  if (location.pathname === "/settings/environments") settingsSection = 'environments';
  else if (location.pathname === "/settings/data-controls") settingsSection = 'data-controls';
  else if (location.pathname === "/settings/general") settingsSection = 'general';

  // Create a new chat session
  const handleNewChat = async () => {
    if (!selectedEnv) return;
    const res = await fetch(`/api/environment/${selectedEnv}/chat_sessions`, {
      method: 'POST',
      credentials: 'include',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ title: null })
    });
    if (res.ok) {
      const session = await res.json();
      setChatSessions(sessions => [session, ...sessions]);
      setSelectedSessionId(session.session_id);
      navigate(`/chat/session/${session.session_id}`);
    }
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
          {/* Settings sidebar */}
          {isSettings ? (
            <>
              <div style={{ marginTop: '1.7rem', marginLeft: '1.2rem', display: 'flex', flexDirection: 'column', alignItems: 'flex-start' }}>
                <button
                  onClick={() => navigate('/')}
                  className="inline-flex items-center justify-center rounded-full focus:outline-none"
                  style={{ background: 'black', cursor: 'pointer', width: '44px', height: '44px', minWidth: '44px', minHeight: '44px', padding: 0, display: 'flex', alignSelf: 'flex-start', marginBottom: 0 }}
                  title="Go to homepage"
                >
                  <span style={{ fontSize: '1.7rem', fontWeight: 900, color: '#22c55e' }}>T</span>
                </button>
                <span style={{ fontSize: '1.2rem', fontWeight: 600, color: '#fff', letterSpacing: '0.5px', marginLeft: 10, marginTop: 18, marginBottom: 16, display: 'block' }}>Settings</span>
                <nav className="flex flex-col gap-2 w-full">
                  <button
                    className={`settings-link py-1 pl-2 pr-2 rounded transition text-left text-sm ${settingsSection === 'general' ? 'bg-green-700 text-white font-bold' : 'text-gray-200'}`}
                    style={{ background: settingsSection === 'general' ? undefined : 'none', minHeight: '32px', fontSize: '0.95rem', paddingLeft: 14, width: '90%', maxWidth: 210, marginLeft: 2 }}
                    onClick={() => navigate('/settings/general')}
                  >General</button>
                  <button
                    className={`settings-link py-1 pl-2 pr-2 rounded transition text-left text-sm ${settingsSection === 'environments' ? 'bg-green-700 text-white font-bold' : 'text-gray-200'}`}
                    style={{ background: settingsSection === 'environments' ? undefined : 'none', minHeight: '32px', fontSize: '0.95rem', paddingLeft: 14, width: '90%', maxWidth: 210, marginLeft: 2 }}
                    onClick={() => navigate('/settings/environments')}
                  >Environments</button>
                  <button
                    className={`settings-link py-1 pl-2 pr-2 rounded transition text-left text-sm ${settingsSection === 'data-controls' ? 'bg-green-700 text-white font-bold' : 'text-gray-200'}`}
                    style={{ background: settingsSection === 'data-controls' ? undefined : 'none', minHeight: '32px', fontSize: '0.95rem', paddingLeft: 14, width: '90%', maxWidth: 210, marginLeft: 2 }}
                    onClick={() => navigate('/settings/data-controls')}
                  >Data controls</button>
                </nav>
              </div>
            </>
          ) : (
            <>
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
              {/* Environment Files List - moved above settings button */}
              {sidebarOpen && !isSettings && selectedEnv && envFiles.length > 0 && (
                <div style={{ marginLeft: '1.2rem', marginTop: '1.5rem', maxHeight: '180px', overflowY: 'auto' }}>
                  <div className="text-m text-gray-400 mb-2 font-semibold flex items-center"><FaListAlt className="mr-2" />Files in environment</div>
                  <ul className="space-y-1">
                    {envFiles.map(file => (
                      <li key={file.file_id} className="flex items-center justify-between text-sm">
                        <button
                          className="truncate max-w-[140px] text-sm text-gray-400 hover:text-green-400 bg-transparent border-none p-0 m-0 text-left cursor-pointer"
                          style={{ background: 'none' }}
                          onClick={() => { setSidebarPreviewFileId(file.file_id); setSidebarPreviewOpen(true); }}
                          title="Preview file"
                        >
                          {file.filename}
                        </button>
                      </li>
                    ))}
                  </ul>
                </div>
              )}
              {/* Chat Sessions List */}
              {sidebarOpen && !isSettings && selectedEnv && (
                <div style={{ marginLeft: '1.2rem', marginTop: '2rem', maxHeight: '180px', overflowY: 'auto' }}>
                  <div className="flex items-center justify-between mb-2">
                    <span className="text-m text-gray-400 font-semibold flex items-center">
                      <FaComments className="mr-2" />Chats
                    </span>
                    <button
                      onClick={handleNewChat}
                      className="px-3 py-1 pr-2 rounded-full text-[#166534] hover:text-[#22c55e] text-xs -ml-1"
                      style={{ background: 'none', border: 'none', boxShadow: 'none', padding: 0, minWidth: 0, minHeight: 0 }}
                      title="New Chat"
                    >
                      <FaPlus />
                    </button>
                  </div>
                  <ul className="space-y-1">
                    {chatSessions.map(session => (
                      <li key={session.session_id}>
                        <button
                          className={`truncate max-w-[180px] text-sm text-left ${selectedSessionId === session.session_id ? 'text-[#22c55e] font-bold' : 'text-[#4ade80] hover:text-[#22c55e]'}`}
                          style={{ width: '100%', background: 'none', border: 'none', boxShadow: 'none', padding: '2px 4px', minHeight: 0, minWidth: 0, borderRadius: 0, textAlign: 'left' }}
                          onClick={() => navigate(`/chat/session/${session.session_id}`)}
                        >
                          {session.title || `Chat ${session.session_id.slice(0, 6)}`}
                        </button>
                      </li>
                    ))}
                  </ul>
                </div>
              )}
              {/* Settings button at the bottom above logout */}
              <div className="flex flex-col items-center mt-auto mb-2">
                <button
                  onClick={() => navigate('/settings/general')}
                  className="flex items-center hover:text-green-400 text-white px-4 py-2 rounded transition"
                  style={{ background: 'none', border: 'none', cursor: 'pointer' }}
                >
                  <svg className="mr-2" width="22" height="22" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
                    <circle cx="12" cy="12" r="3" />
                    <path d="M19.4 15a1.65 1.65 0 0 0 .33 1.82l.06.06a2 2 0 0 1-2.83 2.83l-.06-.06a1.65 1.65 0 0 0-1.82-.33 1.65 1.65 0 0 0-1 1.51V21a2 2 0 0 1-4 0v-.09a1.65 1.65 0 0 0-1-1.51 1.65 1.65 0 0 0-1.82.33l-.06.06a2 2 0 0 1-2.83-2.83l.06-.06a1.65 1.65 0 0 0 .33-1.82 1.65 1.65 0 0 0-1.51-1H3a2 2 0 0 1 0-4h.09a1.65 1.65 0 0 0 1.51-1 1.65 1.65 0 0 0-.33-1.82l-.06-.06a2 2 0 0 1 2.83-2.83l.06.06a1.65 1.65 0 0 0 1.82.33h.09A1.65 1.65 0 0 0 9 3.09V3a2 2 0 0 1 4 0v.09a1.65 1.65 0 0 0 1 1.51h.09a1.65 1.65 0 0 0 1.82-.33l.06-.06a2 2 0 0 1 2.83 2.83l-.06.06a1.65 1.65 0 0 0-.33 1.82v.09a1.65 1.65 0 0 0 1.51 1H21a2 2 0 0 1 0 4h-.09a1.65 1.65 0 0 0-1.51 1z" />
                  </svg>
                  <span>Settings</span>
                </button>
              </div>
            </>
          )}
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
              value={selectedModel}
              onChange={e => setSelectedModel(e.target.value)}
            >
              <option value="toxindex-rap">ToxIndex RAP</option>
              <option value="toxindex-pathway">ToxIndex Pathway</option>
              <option value="toxindex-vanilla">ToxIndex Vanilla</option>
              <option value="toxindex-4th">ToxIndex 4th</option>
              <option value="toxindex-5th">ToxIndex 5th</option>
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
      {auth && (
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
      )}
      {/* Main Content */}
      <main className="flex-1 flex items-center justify-center h-screen">
        {React.Children.map(children, child => {
          if (
            React.isValidElement(child) &&
            (child.type as any).name === 'Dashboard'
          ) {
            return React.cloneElement(
              child as React.ReactElement<any>,
              { 
                selectedModel, 
                setSelectedModel, 
                selectedEnv, 
                setSelectedEnv, 
                environments, 
                setEnvironments,
                selectedSessionId,
                setSelectedSessionId,
                chatSessions,
                handleNewChat
              }
            );
          }
          if (
            React.isValidElement(child) &&
            (child.type as any).name === 'ChatSession'
          ) {
            return React.cloneElement(
              child as React.ReactElement<any>,
              {
                selectedModel,
                setSelectedModel,
                selectedEnv,
                setSelectedEnv,
                environments,
                setEnvironments,
                selectedSessionId,
                setSelectedSessionId,
                chatSessions,
                handleNewChat
              }
            );
          }
          return child;
        })}
      </main>
      <FilePreviewModal
        fileId={sidebarPreviewFileId}
        isOpen={sidebarPreviewOpen}
        onRequestClose={() => setSidebarPreviewOpen(false)}
      />
    </div>
  );
};

export default Layout; 