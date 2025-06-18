import React, { useEffect, useRef, useState } from 'react';
import { useParams, useNavigate } from 'react-router-dom';
import { FaPlus } from 'react-icons/fa';
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';
import UploadCsvModal from './UploadCsvModal';

interface Message {
  role: string;
  content: string;
}

interface Environment {
  environment_id: string;
  title: string;
}

interface ChatSessionProps {
  selectedModel?: string;
  environments: Environment[];
  refreshEnvFiles?: () => void;
  loadingEnvironments: boolean;
  refetchChatSessions?: () => void;
}

const ChatSession: React.FC<ChatSessionProps> = ({ selectedModel, environments, refreshEnvFiles, loadingEnvironments, refetchChatSessions }) => {
  const { sessionId } = useParams<{ sessionId: string }>();
  const navigate = useNavigate();
  const [messages, setMessages] = useState<Message[]>([]);
  const [input, setInput] = useState('');
  const [loading, setLoading] = useState(false);
  const [selectedEnv, setSelectedEnv] = useState<string>('');
  const [showUploadModal, setShowUploadModal] = useState(false);
  const chatWindowRef = useRef<HTMLDivElement>(null);
  const selectRef = useRef<HTMLSelectElement>(null);
  const spanRef = useRef<HTMLSpanElement>(null);
  const [selectWidth] = useState(120);
  const pollingRef = useRef<number | null>(null);

  // Fetch messages for this session
  useEffect(() => {
    if (!sessionId) return;
    fetch(`/api/chat_sessions/${sessionId}/messages`, { credentials: 'include' })
      .then(res => res.json())
      .then(data => setMessages(data.messages || []));
  }, [sessionId]);

  // Scroll to bottom on new messages
  useEffect(() => {
    if (chatWindowRef.current) {
      chatWindowRef.current.scrollTop = chatWindowRef.current.scrollHeight;
    }
  }, [messages]);

  // Helper to check if last message is from assistant
  const lastMessageIsAssistant = (msgs: Message[]) => msgs.length > 0 && msgs[msgs.length - 1].role === 'assistant';

  // Poll for new messages after user sends a message
  const pollForAssistant = () => {
    if (pollingRef.current !== null) clearInterval(pollingRef.current);
    pollingRef.current = window.setInterval(() => {
      fetch(`/api/chat_sessions/${sessionId}/messages`, { credentials: 'include' })
        .then(res => res.json())
        .then(data => {
          setMessages(data.messages || []);
          if (lastMessageIsAssistant(data.messages || [])) {
            if (pollingRef.current !== null) clearInterval(pollingRef.current);
          }
        });
    }, 2000);
  };

  // Start polling for assistant reply when sessionId changes and last message is from user
  useEffect(() => {
    if (!sessionId) return;
    if (messages.length > 0 && messages[messages.length - 1].role === 'user') {
      pollForAssistant();
    }
    // Clean up polling on session change
    return () => { if (pollingRef.current !== null) clearInterval(pollingRef.current); };
  }, [sessionId, messages.length]);

  // Auto-select first environment and restore from localStorage
  useEffect(() => {
    const stored = localStorage.getItem('selectedEnv');
    if (stored && environments.some(e => e.environment_id === stored)) {
      setSelectedEnv(stored);
    } else if (environments.length > 0) {
      setSelectedEnv(environments[0].environment_id);
    }
  }, [environments]);

  useEffect(() => {
    if (selectedEnv) {
      localStorage.setItem('selectedEnv', selectedEnv);
    }
  }, [selectedEnv]);

  const handleSend = async (e: React.FormEvent) => {
    e.preventDefault();
    if (!input.trim() || !sessionId) return;
    const isFirstUserMessage = messages.filter(m => m.role === 'user').length === 0;
    setMessages(msgs => [...msgs, { role: 'user', content: input }]);
    setLoading(true);
    try {
      await fetch(`/api/chat_sessions/${sessionId}/message`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        credentials: 'include',
        body: JSON.stringify({ prompt: input, environment_id: selectedEnv, model: selectedModel }),
      });
      pollForAssistant();
      // Refetch chat sessions if this was the first user message
      if (isFirstUserMessage && typeof refetchChatSessions === 'function') {
        await refetchChatSessions();
        console.log('DEBUG: Called refetchChatSessions after first user message');
      }
    } catch {
      // handle error
    }
    setInput('');
    setLoading(false);
  };

  return (
    <div className="flex flex-col w-full" style={{ height: '100vh', background: 'linear-gradient(135deg, #101614 0%, #1a2a1a 60%, #1a2320 100%)' }}>
      <div className="flex flex-col flex-1 max-w-3xl mx-auto w-full min-h-0" style={{ minHeight: 0, height: '100%' }}>
        {/* Chat ID display at top of scroll box */}
        <div className="w-full text-center text-gray-300 text-sm font-mono" style={{paddingTop: '1.2rem', paddingBottom: '1.2rem'}}>
          {sessionId ? `Chat ID: ${sessionId}` : ''}
        </div>
        {/* Message area */}
        <div
          ref={chatWindowRef}
          className="flex-1 overflow-y-auto px-4 pt-24 pb-6 min-h-0"
          style={{ scrollBehavior: 'smooth', background: 'transparent' }}
        >
          {messages.length === 0 && (
            <div className="text-gray-400 text-center mt-12">No messages yet. Start the conversation!</div>
          )}
          {messages.map((msg, idx) => (
            <div key={idx} className={`flex mb-4 ${msg.role === 'user' ? 'justify-end' : 'justify-start'}`} style={{ width: '100%' }}>
              {msg.role === 'user' ? (
                <div
                  className="rounded-2xl px-5 py-3 shadow-md max-w-[70%] bg-gradient-to-br from-purple-600 to-green-500 text-white ml-auto"
                  style={{
                    textAlign: 'right',
                    fontSize: '1rem',
                    wordBreak: 'break-word',
                    borderBottomRightRadius: 8,
                    borderBottomLeftRadius: 24,
                    borderTopLeftRadius: 24,
                    borderTopRightRadius: 24,
                    marginLeft: 32,
                    marginRight: 0,
                    boxShadow: '0 2px 8px 0 rgba(34,197,94,0.10)'
                  }}
                >
                  <ReactMarkdown remarkPlugins={[remarkGfm]}>{msg.content}</ReactMarkdown>
                </div>
              ) : (
                <div
                  className="text-gray-100 w-full"
                  style={{
                    textAlign: 'left',
                    fontSize: '1rem',
                    wordBreak: 'break-word',
                    marginLeft: 0,
                    lineHeight: 1.8,
                  }}
                >
                  <ReactMarkdown remarkPlugins={[remarkGfm]}>{msg.content}</ReactMarkdown>
                </div>
              )}
            </div>
          ))}
          {/* Show 'Task is running...' if last message is from user */}
          {messages.length > 0 && messages[messages.length-1]?.role === 'user' && (
            <div className="flex justify-start mb-2">
              <div className="rounded-xl px-4 py-2 bg-gray-700 text-white max-w-[70%] text-sm opacity-80">
                Task is running...
              </div>
            </div>
          )}
        </div>
        {/* Input area styled like Dashboard */}
        <form onSubmit={e => {
          e.preventDefault();
          if (selectedEnv === "__add__" || selectedEnv === "__manage__") return;
          handleSend(e);
        }} className="flex flex-col items-center w-full" style={{ maxWidth: '800px', margin: '0 auto', padding: '1.5rem 0', background: 'transparent' }}>
          <div className="relative w-full" style={{ maxWidth: '800px', width: '100%' }}>
            {/* Hidden span to measure width */}
            <span
              ref={spanRef}
              style={{
                position: "absolute",
                visibility: "hidden",
                whiteSpace: "nowrap",
                fontWeight: "bold",
                fontSize: "0.875rem",
                fontFamily: "inherit",
                padding: "0 16px"
              }}
            >
              {selectedEnv === "__add__"
                ? "+ Add environment"
                : selectedEnv === "__manage__"
                  ? "⚙ Manage environments"
                  : (environments ?? []).find(e => e.environment_id === selectedEnv)
                    ? `env - ${(environments ?? []).find(e => e.environment_id === selectedEnv)?.title}`
                    : ""}
            </span>
            <textarea
              rows={3}
              placeholder="Ask me your toxicology question [ Is green tea nephrotoxic? ]"
              value={input}
              onChange={e => setInput(e.target.value)}
              onKeyDown={e => {
                if (e.key === 'Enter' && !e.shiftKey) {
                  e.preventDefault();
                  if (!loading && selectedEnv !== "__add__" && selectedEnv !== "__manage__") handleSend(e);
                }
              }}
              className="w-full pt-4 pb-4 pr-28 pl-8 text-lg rounded-2xl border border-gray-700 bg-gray-900 bg-opacity-70 text-white resize-none min-h-[80px] shadow-2xl focus:outline-none focus:ring-2 focus:ring-green-400 text-left placeholder:text-left"
              style={{ minHeight: 80, fontFamily: 'inherit', width: '100%', boxShadow: '0 8px 32px 0 rgba(34,197,94,0.10)' }}
              autoFocus
              disabled={selectedEnv === "__add__" || selectedEnv === "__manage__"}
            />
            <div className="absolute right-4 bottom-11 flex items-center space-x-2 z-30">
              <button
                className="w-10 h-10 flex items-center justify-center rounded-full shadow-lg hover:bg-gray-200 focus:outline-none focus:ring-2 focus:ring-green-400"
                style={{ padding: 0, borderRadius: '50%', background: 'rgba(255,255,255,0.7)' }}
                title="Upload CSV file"
                onClick={() => setShowUploadModal(true)}
                type="button"
              >
                <FaPlus className="w-5 h-5 text-black" />
              </button>
              <button
                type="submit"
                className="w-10 h-10 flex items-center justify-center rounded-full shadow-lg hover:bg-gray-200 focus:outline-none focus:ring-2 focus:ring-green-400"
                style={{ padding: 0, borderRadius: '50%', background: 'rgba(255,255,255,0.7)' }}
                title="Submit"
                disabled={selectedEnv === "__add__" || selectedEnv === "__manage__"}
              >
                <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" strokeWidth={2.2} stroke="black" className="w-5 h-5">
                  <path strokeLinecap="round" strokeLinejoin="round" d="M5 10l7-7m0 0l7 7m-7-7v18" />
                </svg>
              </button>
            </div>
            <div className="absolute left-4 z-20 flex items-end" style={{ position: 'relative', width: '210px', overflow: 'visible', bottom: '3rem' }}>
              <div className="relative group" style={{ width: '100%', minWidth: '180px', maxWidth: '210px' }}>
                <div className="flex items-center w-full" style={{ position: 'relative', height: '1.7rem', display: 'flex', alignItems: 'center' }}>
                  <select
                    ref={selectRef}
                    className={`font-bold text-white text-sm px-1 py-2 rounded-full appearance-none bg-transparent border-none focus:outline-none transition-all duration-150 group-hover:bg-black group-hover:bg-opacity-60 group-hover:border group-hover:border-gray-700 group-hover:px-4 group-hover:pr-10 group-hover:cursor-pointer focus:bg-black focus:bg-opacity-60 focus:border focus:border-gray-700 focus:px-4 focus:pr-10 w-full ${loadingEnvironments ? 'opacity-50' : ''}`}
                    style={{
                      width: `${selectWidth}px`,
                      minWidth: "50px",
                      maxWidth: "260px",
                      height: '1.7rem',
                      borderRadius: '999px',
                      position: 'relative',
                      zIndex: 1,
                      backgroundColor: 'transparent',
                      cursor: 'pointer',
                      paddingRight: '1.5rem',
                    }}
                    value={selectedEnv}
                    onChange={e => {
                      if (e.target.value === "__add__") {
                        navigate("/settings/environments/create");
                      } else if (e.target.value === "__manage__") {
                        navigate("/settings/environments");
                      } else {
                        setSelectedEnv(e.target.value);
                      }
                    }}
                    disabled={loadingEnvironments}
                  >
                    {(environments ?? []).length > 0 && (environments ?? []).map(env => (
                      <option key={env.environment_id} value={env.environment_id} style={{ paddingLeft: '1rem' }}>
                        {`env - ${env.title}`}
                      </option>
                    ))}
                    <option value="__add__" style={loadingEnvironments ? { opacity: 0.5 } : {}}>
                      {loadingEnvironments ? '+ Add environment' : '+ Add environment'}
                    </option>
                    <option value="__manage__" style={loadingEnvironments ? { opacity: 0.5 } : {}}>&#9881; Manage environments</option>
                  </select>
                </div>
              </div>
            </div>
          </div>
        </form>
      </div>
      <UploadCsvModal
        open={showUploadModal}
        onClose={() => setShowUploadModal(false)}
        environments={environments}
        defaultEnvId={selectedEnv}
        onUploadSuccess={() => {
          setShowUploadModal(false);
          if (refreshEnvFiles) refreshEnvFiles();
        }}
        refreshEnvFiles={refreshEnvFiles}
      />
    </div>
  );
};

export default ChatSession; 