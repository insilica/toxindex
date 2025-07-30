import React, { useEffect, useState, useRef } from "react";
import { useNavigate, useLocation } from "react-router-dom";
import FilePreviewModal from './shared/FilePreviewModal';
import { FaComments, FaPlus, FaListAlt, FaFileCsv, FaFileAlt, FaFileCode, FaDatabase, FaFileImage, FaFile, FaEllipsisH, FaUsers, FaCog, FaServer, FaShieldAlt, FaUser, FaSignOutAlt } from 'react-icons/fa';
import { createPortal } from "react-dom";
import { useEnvironment } from "../context/EnvironmentContext";
import { useChatSession } from "../context/ChatSessionContext";
import { useModel } from "../context/ModelContext";
import { useAdmin } from "../hooks/useAdmin";
import logoUrl from '../assets/logo.svg';

const Layout: React.FC<{ children: React.ReactNode }> = ({ children }) => {
  console.log("Layout mounted");
  const navigate = useNavigate();
  const location = useLocation();
  const [auth, setAuth] = useState<null | boolean>(null);
  const [sidebarOpen, setSidebarOpen] = useState(false);
  const [sidebarClosing, setSidebarClosing] = useState(false);
  const [profileOpen, setProfileOpen] = useState(false);
  const profileRef = useRef<HTMLDivElement>(null);
  const [user, setUser] = useState<{ email?: string; user_id?: string } | null>(null);
  const { selectedModel, setSelectedModel } = useModel();
  const [envFiles, setEnvFiles] = useState<{ file_id: string; filename: string }[]>([]);
  const [sidebarPreviewFileId, setSidebarPreviewFileId] = useState<string | null>(null);
  const [sidebarPreviewOpen, setSidebarPreviewOpen] = useState(false);

  const { refetchChatSessions, selectedSessionId, setSelectedSessionId, chatSessions} = useChatSession();
  const { environments, setEnvironments, selectedEnv, setSelectedEnv } = useEnvironment(); 
  const [openChatMenu, setOpenChatMenu] = useState<string | null>(null);
  const menuButtonRefs = useRef<{ [key: string]: React.RefObject<HTMLButtonElement> }>({});
  const { isAdmin } = useAdmin();

  // Show back arrow on environment and task detail pages
  // const showBackArrow = location.pathname.startsWith('/environments') || location.pathname.startsWith('/settings') || location.pathname.startsWith('/chat/session');
  // const [isBackHover, setIsBackHover] = useState(false);

  useEffect(() => {
    fetch("/api/users/me", { credentials: "include", cache: "no-store" })
      .then(res => {
        if (res.status === 200) {
          setAuth(true);
          return res.json();
        } else {
          setAuth(false);
          return null;
        }
      })
      .then(data => setUser(data))
      .catch(() => {
        setAuth(false);
        setUser(null);
      });
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
      .then(data => setEnvFiles(data.files || []))
      .catch(() => setEnvFiles([]));
  }, [selectedEnv]);

  // Fetch environments and set default selectedEnv (URL param takes priority)
  useEffect(() => {
    fetch("/api/environments", { credentials: "include" })
      .then(res => res.json())
      .then(data => {
        setEnvironments(data.environments || []);
        const params = new URLSearchParams(location.search);
        const envParam = params.get("env");
        if (
          data.environments &&
          data.environments.length > 0 &&
          envParam &&
          data.environments.some((e: { environment_id: string }) => e.environment_id === envParam)
        ) {
          if (selectedEnv !== envParam) setSelectedEnv(envParam);
        } else if (data.environments && data.environments.length > 0) {
          if (selectedEnv !== data.environments[0].environment_id) setSelectedEnv(data.environments[0].environment_id);
        } else {
          if (selectedEnv !== "__add__") setSelectedEnv("__add__");
        }
      })
      .catch(() => {
        setEnvironments([]);
        if (selectedEnv !== "__add__") setSelectedEnv("__add__");
      });
  }, [location.search]);

  // Robust two-way sync between selectedEnv and the URL
  useEffect(() => {
    const params = new URLSearchParams(location.search);
    const envParam = params.get("env");
    const workflowParam = params.get("workflow");

    // If selectedEnv is set and different from URL, update the URL
    if (
      selectedEnv &&
      envParam !== selectedEnv
    ) {
      params.set("env", selectedEnv);
      navigate({ search: params.toString() }, { replace: true });
    }

    // Only set selectedModel from URL if it's not set yet (default is 'toxindex-rap')
    if (
      workflowParam &&
      workflowParam !== selectedModel &&
      (selectedModel === undefined )
    ) {
      setSelectedModel(workflowParam);
      return;
    }

    // If selectedModel is set and different from URL, update the URL
    if (
      selectedModel &&
      workflowParam !== selectedModel
    ) {
      params.set("workflow", selectedModel);
      navigate({ search: params.toString() }, { replace: true });
    }
  }, [environments, selectedEnv, selectedModel, location.search, navigate, setSelectedEnv, setSelectedModel]);

  // Fetch chat sessions (flat, not by environment)
  useEffect(() => {
    refetchChatSessions();
  }, [refetchChatSessions]);

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
  const isAdminPage = location.pathname.startsWith("/admin");
  const isGeneralSettings = location.pathname === "/settings/general";
  let settingsSection: 'general' | 'environments' | 'data-controls' | 'admin' = 'general';
  if (location.pathname === "/settings/environments") settingsSection = 'environments';
  else if (location.pathname === "/settings/data-controls") settingsSection = 'data-controls';
  else if (location.pathname === "/settings/general") settingsSection = 'general';
  else if (location.pathname === "/admin/users") settingsSection = 'admin';

  // Open sidebar by default on general settings page
  useEffect(() => {
    if (isGeneralSettings && !sidebarOpen && !sidebarClosing) {
      setSidebarOpen(true);
    }
  }, [isGeneralSettings]); // Remove sidebarOpen and sidebarClosing from dependencies

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
    setSelectedSessionId(null);
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
          className={`bg-black shadow-sm flex flex-col justify-between min-h-screen relative transition-all duration-500 ease-in-out ${(sidebarOpen || sidebarClosing) ? 'w-64 opacity-100' : 'w-0 opacity-0 pointer-events-none'}`}
          style={{ overflow: 'hidden', position: 'relative', zIndex: 40, padding: 0, margin: 0, borderRight: 'none' }}
        >
          {/* Close button */}
          {(sidebarOpen || sidebarClosing) && (
            <button
              className={`absolute transition-all duration-500 w-9 h-9 flex items-center justify-center rounded-full text-gray-300 hover:text-green-400 shadow`}
              style={{
                background: 'none',
                top: '2rem',
                right: '1.5rem',
                zIndex: 50,
                minWidth: '36px',
                minHeight: '36px',
                padding: 0,
                opacity: sidebarOpen || sidebarClosing ? 1 : 0,
                pointerEvents: sidebarOpen || sidebarClosing ? 'auto' : 'none',
                transform: sidebarOpen && !sidebarClosing ? 'translateX(0)' : 'translateX(-12rem)',
                transition: 'all 0.5s cubic-bezier(0.4,0,0.2,1)'
              }}
              onClick={() => {
                setSidebarClosing(true);
                setSidebarOpen(false);
                setTimeout(() => {
                  setSidebarClosing(false);
                }, 500);
              }}
              title="Collapse sidebar"
            >
              <svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5" strokeLinecap="round" strokeLinejoin="round">
                <polyline points="15 18 9 12 15 6" />
              </svg>
            </button>
          )}
          {/* Settings sidebar */}
          {(isSettings || isAdminPage) ? (
            <>
              <div style={{ marginTop: '1.7rem', marginLeft: '1.2rem', display: 'flex', flexDirection: 'column', alignItems: 'flex-start' }}>
                <button
                  onClick={() => navigate('/')}
                  className="inline-flex items-center justify-center rounded-full focus:outline-none"
                  style={{ background: 'none', cursor: 'pointer', width: '44px', height: '44px', minWidth: '44px', minHeight: '44px', padding: 0, display: 'flex', alignSelf: 'flex-start', marginBottom: 0 }}
                  title="Go to homepage"
                >
                  <img src={logoUrl} alt="Toxindex logo" width={32} height={32} style={{ display: 'block' }} />
                </button>
                <span style={{ fontSize: '1.2rem', fontWeight: 600, color: '#fff', letterSpacing: '0.5px', marginLeft: 10, marginTop: 18, marginBottom: 16, display: 'block' }}>
                  Settings
                </span>
                <nav className="flex flex-col gap-2 w-full">
                  <button
                    className={`settings-link py-1 pl-2 pr-2 rounded transition text-left text-sm ${settingsSection === 'general' ? '!bg-green-700 text-white font-bold' : 'text-gray-200'}`}
                    style={{ background: settingsSection === 'general' ? undefined : 'none', minHeight: '32px', fontSize: '0.95rem', paddingLeft: 14, width: '90%', maxWidth: 210, marginLeft: 2 }}
                    onClick={() => navigate('/settings/general')}
                  >
                    <FaCog className="inline mr-2" />
                    General
                  </button>
                  <button
                    className={`settings-link py-1 pl-2 pr-2 rounded transition text-left text-sm ${settingsSection === 'environments' ? '!bg-green-700 text-white font-bold' : 'text-gray-200'}`}
                    style={{ background: settingsSection === 'environments' ? undefined : 'none', minHeight: '32px', fontSize: '0.95rem', paddingLeft: 14, width: '90%', maxWidth: 210, marginLeft: 2 }}
                    onClick={() => navigate('/settings/environments')}
                  >
                    <FaServer className="inline mr-2" />
                    Environments
                  </button>
                  <button
                    className={`settings-link py-1 pl-2 pr-2 rounded transition text-left text-sm ${settingsSection === 'data-controls' ? '!bg-green-700 text-white font-bold' : 'text-gray-200'}`}
                    style={{ background: settingsSection === 'data-controls' ? undefined : 'none', minHeight: '32px', fontSize: '0.95rem', paddingLeft: 14, width: '90%', maxWidth: 210, marginLeft: 2 }}
                    onClick={() => navigate('/settings/data-controls')}
                  >
                    <FaShieldAlt className="inline mr-2" />
                    Data controls
                  </button>
                  {isAdmin && (
                    <button
                      className={`settings-link py-1 pl-2 pr-2 rounded transition text-left text-sm ${settingsSection === 'admin' ? '!bg-green-700 text-white font-bold' : 'text-gray-200'}`}
                      style={{ background: settingsSection === 'admin' ? undefined : 'none', minHeight: '32px', fontSize: '0.95rem', paddingLeft: 14, width: '90%', maxWidth: 210, marginLeft: 2 }}
                      onClick={() => navigate('/admin/users')}
                    >
                      <FaUsers className="inline mr-2" />
                      Admin
                    </button>
                  )}
                </nav>
              </div>
            </>
          ) : (
            <>
              <div className="flex items-center mb-6" style={{ marginTop: '1.7rem', marginLeft: '1.2rem' }}>
                <button
                  onClick={() => navigate('/')}
                  className="inline-flex items-center justify-center rounded-full focus:outline-none"
                  style={{ background: 'none', cursor: 'pointer', width: '44px', height: '44px', minWidth: '44px', minHeight: '44px', padding: 0, display: 'flex', alignSelf: 'flex-start' }}
                  title="Go to homepage"
                >
                  <img src={logoUrl} alt="Toxindex logo" width={32} height={32} style={{ display: 'block' }} />
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
                      onClick={handleNewChat}
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
                      {chatSessions.map((session: any) => (
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
            className="fixed z-50 w-9 h-9 flex items-center justify-center rounded-full text-gray-300 hover:text-green-400 shadow transition group"
            style={{ top: '1.95rem', left: '1.45rem', minWidth: '32px', minHeight: '32px', position: 'fixed', background: 'none', padding: 0 }}
            onClick={() => setSidebarOpen(true)}
            title="Open sidebar"
          >
            <span className="absolute w-9 h-9 flex items-center justify-center transition-opacity duration-200 opacity-100 group-hover:opacity-0">
              <img src={logoUrl} alt="Toxindex logo" width={32} height={32} style={{ display: 'block' }} />
            </span>
            <span className="absolute w-9 h-9 flex items-center justify-center transition-opacity duration-200 opacity-0 group-hover:opacity-100">
              <svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5" strokeLinecap="round" strokeLinejoin="round">
                <polyline points="9 18 15 12 9 6" />
              </svg>
            </span>
          </button>
        )}
      </div>
      {/* User profile button at top right, aligned with others */}
      {auth && (
        <div style={{ position: 'absolute', top: '2rem', right: '2.5rem', zIndex: 40 }} ref={profileRef}>
          <button
            onClick={() => setProfileOpen(v => !v)}
            className="w-11 h-11 flex items-center justify-center rounded-full !bg-gray-900 !bg-opacity-10 text-white hover:!bg-gray-800 focus:outline-none shadow"
            style={{ padding: 0 }}
            title="User profile"
          >
            <svg width="28" height="28" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth="2" className="text-gray-300">
              <circle cx="12" cy="8" r="4" />
              <path d="M4 20c0-2.5 3.5-4 8-4s8 1.5 8 4" />
            </svg>
          </button>
          {profileOpen && (
            <div className="absolute right-0 mt-2 w-56 bg-neutral-200 rounded-lg shadow-lg py-2 text-gray-900 border border-gray-200">
              {user?.user_id ? (
                <button
                  onClick={() => {
                    navigate(`/user/${user.user_id}`);
                    setProfileOpen(false);
                  }}
                  className="font-medium text-neutral-900 hover:text-green-500 transition-colors cursor-pointer w-full text-left px-4 py-2 border-b border-gray-100 text-sm"
                  style={{ background: 'none', border: 'none' }}
                >
                  <FaUser className="inline mr-2" />
                  {user.email || "Logged in"}
                </button>
              ) : (
                <div className="font-medium flex items-center px-4 py-2 border-b border-gray-100 text-sm">
                  <FaUser className="inline mr-2" />
                  {user?.email || "Logged in"}
                </div>
              )}
              <button
                onClick={() => {
                  navigate('/settings/general');
                  setProfileOpen(false);
                }}
                className="w-full text-left hover:text-purple-500 px-4 py-2 text-gray-900 text-sm"
                style={{ background: 'none', border: 'none' }}
              >
                <FaCog className="inline mr-2" />
                Settings
              </button>
              <button
                onClick={handleLogout}
                className="w-full text-left hover:text-red-500 px-4 py-2 text-gray-900 text-sm"
                style={{ background: 'none', border: 'none' }}
              >
                <FaSignOutAlt className="inline mr-2" />
                Logout
              </button>
            </div>
          )}
        </div>
      )}
      
      {/* Main Content */}
      <main className="flex-1 flex items-center justify-center h-screen !bg-gray-900 relative">
        {/* {children} */}
        {(
          React.isValidElement(children) &&
          ((children.type as any).displayName === "EnvironmentDetails" || (children.type as any).name === "EnvironmentDetails")
        )
          ? React.cloneElement(children as React.ReactElement<any>, { paddingClass: sidebarOpen ? "px-90" : "px-125" })
          : children}
        

      </main>
      <FilePreviewModal
        fileId={sidebarPreviewFileId}
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