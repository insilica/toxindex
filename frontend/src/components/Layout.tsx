import React, { useEffect, useState, useRef } from "react";
import { useNavigate, useLocation } from "react-router-dom";
import FilePreviewModal from './FilePreviewModal';
import { FaComments, FaPlus, FaListAlt, FaFileCsv, FaFileAlt, FaFileCode, FaDatabase, FaFileImage, FaFile, FaEllipsisH } from 'react-icons/fa';
import { createPortal } from "react-dom";

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
  const [openChatMenu, setOpenChatMenu] = useState<string | null>(null);
  const menuButtonRefs = useRef<{ [key: string]: React.RefObject<HTMLButtonElement> }>({});

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
    if (!selectedEnv) {
      setEnvFiles([]); // Clear files if no environment selected
      return;
    }
    fetch(`/api/environments/${selectedEnv}/files`, { credentials: 'include' })
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

  // Fetch chat sessions (flat, not by environment)
  useEffect(() => {
    fetch(`/api/chat_sessions`, { credentials: 'include' })
      .then(res => res.json())
      .then(data => {
        setChatSessions(data.sessions || []);
        if (data.sessions && data.sessions.length > 0) {
          setSelectedSessionId(data.sessions[0].session_id);
        } else {
          setSelectedSessionId(null);
        }
      });
  }, []);

  const handleLogout = () => {
    fetch("/api/auth/logout", {
      method: "POST",
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

  // Helper to refetch chat sessions (flat)
  const refetchChatSessions = async () => {
    const res = await fetch(`/api/chat_sessions`, { credentials: 'include' });
    const data = await res.json();
    setChatSessions(data.sessions || []);
    console.log('DEBUG: Sidebar chat sessions updated', data.sessions);
  };

  // Create a new chat session (flat)
  const handleNewChat = async () => {
    const res = await fetch(`/api/chat_sessions`, {
      method: 'POST',
      credentials: 'include',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ title: null })
    });
    if (res.ok) {
      const session = await res.json();
      await refetchChatSessions();
      setSelectedSessionId(session.session_id);
      navigate(`/chat/session/${session.session_id}`);
    }
  };

  // Add a function to refresh the file list for the current environment
  const refreshEnvFiles = () => {
    if (!selectedEnv) {
      setEnvFiles([]);
      return;
    }
    fetch(`/api/environments/${selectedEnv}/files`, { credentials: 'include' })
      .then(res => res.json())
      .then(data => setEnvFiles(data.files || []));
  };

  function getFileIcon(filename: string) {
    const ext = filename.split('.').pop()?.toLowerCase();
    if (!ext) return <FaFile />;
    if (ext === 'csv') return <FaFileCsv className="text-green-400" title="CSV file" />;
    if (ext === 'txt') return <FaFileAlt className="text-gray-400" title="Text file" />;
    if (ext === 'json') return <FaFileCode className="text-yellow-400" title="JSON file" />;
    if (ext === 'parquet') return <FaDatabase className="text-blue-400" title="Parquet file" />;
    if (["png","jpg","jpeg","gif","bmp","webp"].includes(ext)) return <FaFileImage className="text-purple-400" title="Image file" />;
    return <FaFile className="text-gray-500" title="File" />;
  }

  // Handler to delete a chat session
  const handleDeleteChat = async (sessionId: string) => {
    if (!window.confirm('Are you sure you want to delete this chat?')) return;
    await fetch(`/api/chat_sessions/${sessionId}`, {
      method: 'DELETE',
      credentials: 'include',
    });
    await refetchChatSessions();
    setOpenChatMenu(null);
    // If the deleted chat was selected, clear selection
    if (selectedSessionId === sessionId) {
      setSelectedSessionId(null);
    }
  };

  // Handler to rename a chat session
  const handleRenameChat = async (sessionId: string, currentTitle: string) => {
    const newTitle = window.prompt('Enter new chat name:', currentTitle);
    if (!newTitle || newTitle.trim() === '' || newTitle === currentTitle) return;
    await fetch(`/api/chat_sessions/${sessionId}`, {
      method: 'PATCH',
      credentials: 'include',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ title: newTitle })
    });
    await refetchChatSessions();
    setOpenChatMenu(null);
  };

  // Add a click-away listener to close the dropdown
  useEffect(() => {
    function handleClick(e: MouseEvent) {
      if (openChatMenu && !(e.target as HTMLElement).closest('.group')) {
        setOpenChatMenu(null);
      }
    }
    if (openChatMenu) document.addEventListener('mousedown', handleClick);
    return () => document.removeEventListener('mousedown', handleClick);
  }, [openChatMenu]);

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
              className={`absolute transition-all duration-300 w-9 h-9 flex items-center justify-center rounded-full text-gray-300 hover:text-green-400 shadow`}
              style={{
                background: 'none',
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
              <svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5" strokeLinecap="round" strokeLinejoin="round">
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
              {/* Environment Files List - always show if environment is selected */}
              {sidebarOpen && !isSettings && selectedEnv && (
                <div style={{ marginLeft: '1.2rem', marginTop: '1.5rem', maxHeight: '180px', overflowY: 'auto' }}>
                  <div className="text-m text-gray-400 mb-2 font-semibold flex items-center"><FaListAlt className="mr-2" />Files in environment</div>
                  {envFiles.length === 0 ? (
                    <div className="text-gray-500 text-sm">No files yet.</div>
                  ) : (
                    <ul className="space-y-1">
                      {envFiles.map(file => (
                        <li key={file.file_id} className="flex items-center justify-between text-sm">
                          <button
                            className="truncate max-w-[250px] text-sm text-gray-400 hover:text-green-400 bg-transparent border-none p-0 m-0 text-left cursor-pointer flex items-center gap-1"
                            style={{ background: 'none' }}
                            onClick={() => { setSidebarPreviewFileId(file.file_id); setSidebarPreviewOpen(true); }}
                            title="Preview file"
                          >
                            {getFileIcon(file.filename)}
                            {file.filename}
                          </button>
                        </li>
                      ))}
                    </ul>
                  )}
                </div>
              )}
              {/* Chat Sessions List */}
              {sidebarOpen && !isSettings && selectedEnv && (
                <div style={{ marginLeft: '1.2rem', marginTop: '2rem', maxHeight: '600px', overflowY: 'auto' }}>
                  <div className="flex items-center justify-between mb-2">
                    <span className="text-m text-gray-400 font-semibold flex items-center">
                      <FaComments className="mr-2" />Chats
                    </span>
                    <button
                      onClick={() => {
                        console.log('handleNewChat clicked, selectedEnv:', selectedEnv);
                        handleNewChat();
                      }}
                      disabled={!selectedEnv || selectedEnv === "__add__" || selectedEnv === "__manage__"}
                      className="px-3 py-1 pr-2 rounded-full text-[#166534] hover:text-[#22c55e] text-xs -ml-1"
                      style={{ background: 'none', border: 'none', boxShadow: 'none', padding: 0, minWidth: 0, minHeight: 0, paddingRight: 35}}
                      title="New Chat"
                    >
                      <FaPlus />
                    </button>
                  </div>
                  {!selectedEnv && <div style={{color: 'red', fontSize: '0.8rem'}}>No environment selected</div>}
                  {chatSessions.length === 0 ? (
                    <div className="text-gray-500 text-sm">No chats yet.</div>
                  ) : (
                    <ul className="space-y-1">
                      {chatSessions.map(session => (
                        <li key={session.session_id} className="relative group">
                          <button
                            className={`truncate max-w-[180px] text-sm text-left ${selectedSessionId === session.session_id ? 'text-[#22c55e] font-bold' : 'text-[#4ade80] hover:text-[#22c55e]'}`}
                            style={{ width: 'calc(100% - 32px)', background: 'none', border: 'none', boxShadow: 'none', padding: '2px 4px', minHeight: 0, minWidth: 0, borderRadius: 0, textAlign: 'left', display: 'inline-block' }}
                            onClick={() => navigate(`/chat/session/${session.session_id}`)}
                          >
                            {(session.title && session.title.length > 0)
                              ? session.title.slice(0, 40)
                              : `Chat ${session.session_id.slice(0, 6)}`}
                          </button>
                          <button
                            ref={(() => {
                              if (!menuButtonRefs.current[session.session_id]) {
                                menuButtonRefs.current[session.session_id] = React.createRef<HTMLButtonElement>() as React.RefObject<HTMLButtonElement>;
                              }
                              return menuButtonRefs.current[session.session_id] as React.RefObject<HTMLButtonElement>;
                            })()}
                            className="ml-1 p-1 text-gray-400 hover:text-green-400 focus:outline-none"
                            style={{
                              background: 'none',
                              boxShadow: 'none',
                              border: 'none',
                              position: 'absolute',
                              right: 0,
                              top: 0,
                              zIndex: 10,
                              padding: 2,
                              paddingRight: 25,
                              minHeight: 0,
                              borderRadius: 4,
                              transition: 'color 0.15s'
                            }}
                            onClick={e => {
                              e.stopPropagation();
                              setOpenChatMenu(openChatMenu === session.session_id ? null : session.session_id);
                            }}
                            title="Chat options"
                          >
                            <FaEllipsisH size={12} />
                          </button>
                          {openChatMenu === session.session_id && menuButtonRefs.current[session.session_id] && <ChatMenuPortal anchorRef={menuButtonRefs.current[session.session_id] as React.RefObject<HTMLButtonElement | null>} onRename={() => handleRenameChat(session.session_id, session.title || '')} onDelete={() => handleDeleteChat(session.session_id)} onClose={() => setOpenChatMenu(null)} />}
                        </li>
                      ))}
                    </ul>
                  )}
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
            className="fixed left-2 z-50 w-9 h-9 flex items-center justify-center rounded-full text-gray-300 hover:text-green-400 shadow transition"
            style={{ top: '2rem', minWidth: '36px', minHeight: '36px', padding: 0, position: 'fixed', background: 'none', paddingLeft: 10}}
            onClick={() => setSidebarOpen(true)}
            title="Open sidebar"
          >
            <svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5" strokeLinecap="round" strokeLinejoin="round">
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
                handleNewChat,
                refetchChatSessions
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
                handleNewChat,
                refetchChatSessions,
                refreshEnvFiles
              }
            );
          }
          // Pass refreshEnvFiles to EnvironmentDetails
          if (
            React.isValidElement(child) &&
            (child.type as any).name === 'EnvironmentDetails'
          ) {
            return React.cloneElement(
              child as React.ReactElement<any>,
              {
                refreshEnvFiles
              }
            );
          }
          return child;
        })}
      </main>
      <FilePreviewModal
        fileId={sidebarPreviewFileId}
        envId={selectedEnv}
        isOpen={sidebarPreviewOpen}
        onRequestClose={() => setSidebarPreviewOpen(false)}
      />
    </div>
  );
};

const ChatMenuPortal: React.FC<{ anchorRef: React.RefObject<HTMLButtonElement | null>, onRename: () => void, onDelete: () => void, onClose: () => void }> = ({ anchorRef, onRename, onDelete, onClose }) => {
  const [pos, setPos] = useState<{ top: number, left: number }>({ top: 0, left: 0 });
  const menuRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    function updatePos() {
      if (anchorRef && anchorRef.current) {
        const rect = anchorRef.current.getBoundingClientRect();
        setPos({
          top: rect.bottom + window.scrollY,
          left: rect.left + window.scrollX
        });
      }
    }
    updatePos();
    window.addEventListener('resize', updatePos);
    window.addEventListener('scroll', updatePos, true);
    return () => {
      window.removeEventListener('resize', updatePos);
      window.removeEventListener('scroll', updatePos, true);
    };
  }, [anchorRef]);

  useEffect(() => {
    function handleClick(e: MouseEvent) {
      if (
        menuRef.current &&
        !menuRef.current.contains(e.target as Node) &&
        anchorRef.current &&
        !anchorRef.current.contains(e.target as Node)
      ) {
        console.log('onClose from document click');
        onClose();
      }
    }
    document.addEventListener('mousedown', handleClick);
    return () => document.removeEventListener('mousedown', handleClick);
  }, [onClose, anchorRef]);

  return createPortal(
    <div ref={menuRef} className="bg-white rounded shadow-lg z-50 border border-gray-200" style={{ minWidth: 120, position: 'absolute', top: pos.top, left: pos.left, width: 128 }}>
      <button
        className="block w-full text-left px-4 py-2 text-sm text-gray-700 hover:bg-gray-100"
        onMouseDown={e => {
          e.preventDefault();
          e.stopPropagation();
          console.log('Rename clicked');
          onRename();
          console.log('onClose after rename');
          onClose();
        }}
      >
        Rename
      </button>
      <button
        className="block w-full text-left px-4 py-2 text-sm text-red-600 hover:bg-red-50"
        onMouseDown={e => {
          e.preventDefault();
          e.stopPropagation();
          console.log('Delete clicked');
          onDelete();
          console.log('onClose after delete');
          onClose();
        }}
      >
        Delete
      </button>
    </div>,
    document.body
  );
};

export default Layout; 