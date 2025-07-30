import React, { useEffect, useState, useRef } from 'react';
import { useParams } from 'react-router-dom';
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';
import { useEnvironment } from "../context/EnvironmentContext";
import { useModel } from "../context/ModelContext";
import ChatInputBar from './shared/ChatInputBar';
import HomeButton from './shared/HomeButton';
import LoadingSpinner from './shared/LoadingSpinner';
import { getWorkflowId } from './shared/workflows';
import { useSocket } from '../context/SocketContext';

let mountCount = 0;

const ChatSession = () => {
  const { sessionId } = useParams<{ sessionId: string }>();
  const [messages, setMessages] = useState<any[]>([]);
  const [loading, setLoading] = useState(true);
  const [sending, setSending] = useState(false);
  const [input, setInput] = useState('');
  const [error, setError] = useState<string | null>(null);
  const messagesEndRef = useRef<HTMLDivElement>(null);
  const { selectedEnv } = useEnvironment();
  const { selectedModel } = useModel();
  const [fileId, setFileId] = useState<string | undefined>(undefined);
  const [fileName, setFileName] = useState<string | undefined>(undefined);
  const { socket, isConnected } = useSocket();

  console.log("ChatSession mounted", sessionId);

  useEffect(() => {
    mountCount += 1;
    console.log("ChatSession useEffect mount count:", mountCount);
  }, []);

  useEffect(() => {
    if (!socket || !sessionId) return;

    // Only join room if socket is connected
    if (isConnected) {
      // Join the room
      console.log('[SocketIO] Emitting join_chat_session', sessionId);
      socket.emit('join_chat_session', { session_id: sessionId });
    }

    // Handler for new_message
    const handler = (msg: any) => {
      console.log('[SocketIO] Received new_message:', msg);
      setMessages(prev => [...prev, msg]);
    };
    socket.on('new_message', handler);

    // Cleanup handler only (not the socket)
    return () => {
      socket.off('new_message', handler);
    };
  }, [sessionId, socket, isConnected]);

  useEffect(() => {
    if (!sessionId) return;
    setLoading(true);
    fetch(`/api/chat_sessions/${sessionId}/messages`, { credentials: 'include' })
      .then(res => res.json())
      .then(data => {
        setMessages(data.messages || []);
      })
      .catch(err => {
        setMessages([]);
        setError("Failed to load messages.");
        console.error("Error loading messages:", err);
      })
      .finally(() => setLoading(false));
  }, [sessionId]);

  const scrollToBottom = () => {
    messagesEndRef.current?.scrollIntoView({ behavior: "smooth" });
  };

  useEffect(() => {
    scrollToBottom();
  }, [messages]);

  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    if (!input.trim() || sending) return;
    const workflow_id = getWorkflowId(selectedModel);
    // Block submit if workflow 4 and no file selected
    if (workflow_id === 4 && !fileId) {
      setError("You must select a file for this workflow.");
      return;
    }
    setSending(true);
    setError(null);
    if (workflow_id === 0) {
      setError("This workflow is not yet supported.");
      setSending(false);
      return;
    }
    try {
      const res = await fetch('/api/tasks', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          message: input.trim(),
          workflow: workflow_id,
          environment_id: selectedEnv || undefined,
          sid: sessionId,
          file_id: fileId,
        }),
        credentials: 'include'
      });
      const data = await res.json();
      if (data.error) {
        setError(data.error);
      } else {
        setMessages(prev => [...prev, { content: input.trim(), role: 'user' }]);
        setInput('');
        setFileId(undefined);
        setFileName(undefined);
      }
    } catch (error) {
      setError('Error sending message.');
      console.error('Error sending message:', error);
    } finally {
      setSending(false);
    }
  };

  return (
    <div className="min-h-screen flex flex-col flex-1 relative" style={{ background: 'linear-gradient(135deg, #101614 0%, #1a2a1a 60%, #1a2320 100%)' }}>
      <div className="flex-1 flex flex-col max-w-5xl mx-auto w-full px-4 pt-20 py-0">
        {/* Chat history box */}
        <div
          className="flex-1 overflow-y-auto w-full max-w-5xl mx-auto px-20 py-0 bg-transparent rounded-lg shadow-inner"
          style={{
            scrollBehavior: 'smooth',
            background: 'transparent',
            minHeight: 0, // allow flex child to shrink
            maxHeight: '75vh', // limit height so only this box scrolls
          }}
        >
          {loading ? (
            <LoadingSpinner text="Loading messages..." />
          ) : (
            <>
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
                      <ReactMarkdown
                        remarkPlugins={[remarkGfm]}
                        components={{
                          td: ({node, ...props}) => <td className="markdown-table-cell" {...props} />
                        }}
                      >
                        {msg.content}
                      </ReactMarkdown>
                    </div>
                  )}
                </div>
              ))}
              {messages.length > 0 && messages[messages.length-1]?.role === 'user' && (
                <div className="flex justify-start mb-2">
                  <div className="rounded-xl px-4 py-2 bg-gray-700 text-white max-w-[70%] text-sm opacity-80">
                    Task is running...
                  </div>
                </div>
              )}
              <div ref={messagesEndRef} />
            </>
          )}
        </div>
        {/* Selected file display above input bar, left-aligned */}
        <div className="w-full flex flex-col items-center pb-8 relative">
          <div style={{ width: '100%', maxWidth: 800 }}>
            {fileId && (
              <div className="bg-gray-900/10 text-green-200 px-3 py-0 rounded-full text-sm font-medium flex items-center gap-1 mb-2" style={{ minHeight: 40 }}>
                <span className="font-mono font-semibold text-base z-13">Selected file:</span>
                <span className="font-mono truncate max-w-xs text-base z-13">{fileName ? (fileName.length > 24 ? fileName.slice(0, 20) + '...' : fileName) : fileId.slice(0, 10) + '...'}</span>
                <button onClick={() => { setFileId(undefined); setFileName(undefined); }} className="ml-0 text-red-300 hover:text-red-500 z-13" title="Clear file">×</button>
              </div>
            )}
            <ChatInputBar
              value={input}
              onChange={setInput}
              onSubmit={handleSubmit}
              uploading={sending}
              error={error}
              onFilePick={(id: string, name?: string) => {
                setFileId(id);
                setFileName(name);
              }}
            />
          </div>
          
          {/* Home button positioned relative to page */}
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
      </div>
    </div>
  );
};

export default ChatSession; 